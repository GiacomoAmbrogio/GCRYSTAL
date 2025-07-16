Module interface_MPI
Use base_numbers
Use MPI

!===========================================================
!    @author: Giacomo Ambrogio       date: Jun 2024
!
!    DASH - Device Accelerated Solutions for HPC
!
!    Interface to MPI
!
!    Simplified version for wrapping MPI procedures,
!    nvfortran seesm to have some issues with some of the
!    featurs of MPI_f08, so stick to the 'old' MPI 90
!===========================================================

  Implicit None
  Private

  Integer, Parameter, Public :: NOT_INITIALIZED = -99
  Integer, Parameter, Public :: INITIALIZED     =  1

  Integer, Parameter, Public :: DASH_TYPE_INT  = MPI_INTEGER
  Integer, Parameter, Public :: DASH_TYPE_REAL = MPI_DOUBLE_PRECISION
  Integer, Parameter, Public :: DASH_TYPE_COMP = MPI_DOUBLE_COMPLEX
  Integer, Parameter, Public :: DASH_TYPE_LOG  = MPI_LOGICAL

  Type, Public :: proc_info
    Integer, Private :: comm                                    !! MPI Communicator
    Integer, Private :: size                                    !! Communicator size
    Integer, Private :: rank                                    !! Communicator rank
    Integer, Private :: status   = NOT_INITIALIZED              !! Current object status
    Integer, Private :: my_point = NOT_INITIALIZED              !! Indeces pointer
    Contains
    Generic  ,                  Public  :: sum            => proc_info_sum, proc_info_sum_r1
    Generic  ,                  Public  :: gather         => proc_info_gather,proc_info_gather_r1,proc_info_gather_real
    Procedure,                  Public  :: create         => proc_info_create_from_comm
    Procedure,                  Public  :: get            => proc_info_get
    Procedure,                  Public  :: is_initialized => proc_info_is_initialized
    Procedure, Non_overridable, Public  :: sync           => proc_info_sync
    Procedure, Non_overridable, Public  :: max            => proc_info_max
    Procedure, Non_overridable, Public  :: broadcast      => proc_info_broadcast_r1
    Procedure, Non_overridable, Public  :: my_point_reset => proc_info_mp_reset
    Procedure, Non_overridable, Public  :: my_point_nexy  => proc_info_mp_next
    Procedure, Non_overridable, Private :: proc_info_sum
    Procedure, Non_overridable, Private :: proc_info_sum_r1
    Procedure, Non_overridable, Private :: proc_info_broadcast_r1
    Procedure, Non_overridable, Private :: proc_info_gather
    Procedure, Non_overridable, Private :: proc_info_gather_r1
    Procedure, Non_overridable, Private :: proc_info_gather_real
  End Type

  Type, Public :: data_file_metadata
    Integer      , Private :: record                            !! index of the ks
    Integer(fp64), Private :: record_lenght = INVALID           !! total disk size for the ks (can be larger than needed)
    Integer      , Private :: record_size                       !! how many ks points in this process
    Integer      , Private :: record_type   = POINT_NOT_EXIST   !! type of this ks
    Contains
    Procedure, Public :: create => data_file_metadata_create
  End Type

  Type, Public :: data_file
    Integer                 ,                            Private :: unit
    Type(data_file_metadata), Dimension(:), Allocatable, Private :: metadata
    Integer                 ,                            Private :: status   = NOT_INITIALIZED
    Contains
    Procedure, Public :: init               => data_file_open_init
    Procedure, Public :: initialized        => data_file_initialized
    Procedure, Public :: close              => data_file_close
    Generic  , Public :: data_to_disk       => real_write_to_disk, comp_write_to_disk
    Generic  , Public :: data_from_disk     => real_read_from_disk, comp_read_from_disk
    Procedure, Public :: vect_from_disk     => data_file_vect_real_read_from_disk
    Procedure, Public :: transfer_ownership => data_file_transfer_ownership
    ! Private Implementations
    Procedure, Private :: real_write_to_disk  => data_file_real_write_to_disk
    Procedure, Private :: comp_write_to_disk  => data_file_comp_write_to_disk
    Procedure, Private :: real_read_from_disk => data_file_real_read_from_disk
    procedure, Private :: comp_read_from_disk => data_file_comp_read_from_disk
  End Type

  !! Global MPI communicator (copy of MPI_COMM_WORLD)
  Type(proc_info), Public :: world

  !! --- Utility ---
  Public :: dash_mpi_init_context            !! initialize basic parallelization info
  Public :: dash_mpi_interface_stop_error    !! error check
  !! --- unit/file managment ---
  !Public :: dash_mpi_open_unit
  !Public :: dash_mpi_close_unit
  !! --- k/s distribution ---
  Public :: dash_mpi_assign_batches
  Public :: dash_mpi_assign_batches_manual

  Contains

  !! ******************************************
  !! ============ Public Procedures ===========
  !! ******************************************

    Subroutine dash_mpi_init_context()
    !! Set up the MPI contexts
    !! Creating world communicator
      Call world%create(MPI_COMM_WORLD)
    End Subroutine dash_mpi_init_context

    Subroutine dash_mpi_assign_batches(env,n_batches,my_batch)
    !! Automatic k/s points batches assignment to processes
    !! This procedure should distribute offloding batches to diffrent nodes
    !! equally, and also set device indeces in a round robin fasion to
    !! offloading processes
    Use backend_module, Only: bk_get_visible_devices,bk_set_current_device
      Type(proc_info), Intent(In   ) :: env
      Integer,         Intent(In   ) :: n_batches
      Integer,         Intent(  Out) :: my_batch
      ! Local Variables
      Integer :: node_comm, head_comm
      Integer :: ierr
      Integer :: my_node_rank, my_node_size
      Integer :: node_ID
      Integer :: total_nodes
      Integer :: this_batch
      Integer :: batches_per_node, spear_batches
      Integer :: my_start, my_end
      Integer :: counter
      Integer :: node_devices, my_device_index
      !! create a temporary node communicator
      Call MPI_comm_split_type(env%comm,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,node_comm,ierr)
      Call MPI_comm_rank(node_comm,my_node_rank,ierr)
      Call MPI_comm_size(node_comm,my_node_size,ierr)
      !! Create a communicator for all rank 0 in each node
      Call MPI_comm_split(env%comm,my_node_rank,0,head_comm,ierr)
      !! This gives a unique identifiers for each node, reliable only for my_node_rank = 0
      Call MPI_comm_rank(head_comm,node_ID,ierr)
      If(my_node_rank.eq.0) Then
        total_nodes = 1
      Else
        total_nodes = 0
        node_ID = -9
      Endif
      !! Gathering info on total number of nodes from all processes in env
      Call env%sum(total_nodes)
      !! Getting every process to know its node_ID
      Call MPI_allreduce(MPI_IN_PLACE,node_ID,1,get_type(node_ID),MPI_MAX,node_comm,ierr)
      !! --- Now perform real distribution of batches: ---
      !! How many batches in each node
      batches_per_node = n_batches / total_nodes
      !! How many batches remain
      spear_batches = Mod(n_batches,total_nodes)
      If (node_ID.lt.spear_batches) Then
        my_start = node_ID * (batches_per_node + 1)
        my_end = my_start + batches_per_node
        If ((batches_per_node+1).gt.my_node_size) Then
          Call dash_mpi_interface_stop_error(0,"dash_mpi","too many offloader required")
        Endif
      Else
        my_start = node_ID * batches_per_node + spear_batches
        my_end = my_start + batches_per_node - 1
        If ((batches_per_node).gt.my_node_size) Then
          Call dash_mpi_interface_stop_error(0,"dash_mpi","too many offloader required")
        Endif
      Endif
      !! Assigning batch ID and device ID to processes
      counter = 0
      my_batch = INVALID
      node_devices = bk_get_visible_devices(limit=DEVICE_LIMIT)
      Do this_batch = my_start, my_end
        If(my_node_rank.eq.counter) Then
          my_batch = this_batch
          my_device_index = Mod(counter,node_devices)
          Call bk_set_current_device(my_device_index)
        Endif
        counter = counter + 1
      Enddo
      !! Free MPI communicators
      Call MPI_comm_free(head_comm,ierr)
      Call MPI_comm_free(node_comm,ierr)
    End Subroutine dash_mpi_assign_batches

    Subroutine dash_mpi_assign_batches_manual(env,n_batches,my_batch,rank_list,dev_list)
    Use backend_module, Only: bk_set_current_device
      Type(proc_info),       Intent(In   )           :: env
      Integer,               Intent(In   )           :: n_batches
      Integer,               Intent(  Out)           :: my_batch
      Integer, Dimension(:), Intent(In   )           :: rank_list
      Integer, Dimension(:), Intent(In   ), Optional :: dev_list
      ! Local Variables
      Integer :: i
      Integer :: occurrencies
      Integer :: my_dev_idx
      !! Batches assignment
      my_batch = INVALID
      Loop_on_batches: Do i = 1, n_batches
        occurrencies = Count(rank_list.eq.rank_list(i))
        If(occurrencies.gt.1) Then
          Call dash_mpi_interface_stop_error(0,"dash_mpi","assigning more batches to same process")
        Endif
        If(env%rank.eq.rank_list(i)) Then
          !! should start from zero
          my_batch = i - 1
          !! Devices assignment
          my_dev_idx = 0 !! if not specified all process will use gpu 0
          If(Present(dev_list)) Then
            my_dev_idx = dev_list(i)
          Endif
          Call bk_set_current_device(my_dev_idx,DEVICE_LIMIT)
        Endif
      Enddo Loop_on_batches
    End Subroutine dash_mpi_assign_batches_manual

    Subroutine dash_mpi_interface_stop_error(code,label,message,data)
    !! Print messages based on the error code
    !! Stop execution on code 0 occurence
      Integer         , Intent(In)           :: code
      Character(Len=*), Intent(In)           :: label
      Character(Len=*), Intent(In)           :: message
      Class(*)        , Intent(In), Optional :: data
      ! Local Variables
      Character(Len=26) :: error_data
      Character(Len=40) :: form_err  = '(" ***** ERROR :: ",A," ***** ",A,1X,A)'
      Character(Len=42) :: form_warn = '(" ***** WARNING :: ",A," ***** ",A,1X,A)'
      Character(Len=46) :: form_info = '(" ***** INFORMATION :: ",A," ***** ",A,1X,A)'
      If(Present(data))Then
        Select Type(data)
        Type Is(Integer)
          Write(error_data,'(I0)') data
        Type Is(Logical)
          Write(error_data,'(L1)') data
        Type Is(Real(fp64))
          Write(error_data,'(E18.4)') data
        Type Is(Complex(fp64))
          Write(error_data,'(E12.4,1X,E12.4,"i")') data
        Class Default
          Write(*,*) 'Invalid data type in dash_mpi_interface_stop_error'
          Call dash_stop
        End Select
      Else
        error_data = ''
      Endif
      Select Case(code)
      Case(0)
        Write(6,form_err) label,message,Trim(error_data)
        Call dash_stop
      Case(1)
        If(world%rank.eq.0) Write(*,form_warn) label,message,Trim(error_data)
      Case Default
        If(world%rank.eq.0) Write(*,form_info) label,message,Trim(error_data)
      End Select
    End Subroutine dash_mpi_interface_stop_error

  !! ********************************************
  !! ============ I/O Units Managment ===========
  !! ********************************************
    !Subroutine dash_mpi_open_unit(iunit,name,parallel)
    !!! Used to open units, both parallelized ones and serial ones
    !  Integer         , Intent(In) :: iunit
    !  Character(Len=*), Intent(In) :: name
    !  Logical         , Intent(In) :: parallel
    !  ! Locala Variables
    !  Logical           :: stat
    !  Integer           :: ierr
    !  Logical           :: already_opened
    !  Character(Len=10) :: suffix
    !  stat = world%is_initialized()
    !  Inquire(Unit=iunit,Opened=already_opened)
    !  If(already_opened) Call dash_mpi_interface_stop_error(0,'open_unit_file','opening file already opened')
    !  ierr = 0
    !  If(parallel) Then
    !    Write(suffix,'(A,I0)') '.dash',world%rank
    !    Open(Unit=iunit,File=name//suffix,Iostat=ierr)
    !  Else
    !    If(world%rank.eq.0) Open(Unit=iunit,File=name,Iostat=ierr)
    !  Endif
    !  If(ierr.ne.0) Call dash_mpi_interface_stop_error(0,'open_unit_file','opening unit:',iunit)
    !End Subroutine dash_mpi_open_unit

    !Subroutine dash_mpi_close_unit(iunit)
    !  Integer, Intent(In) :: iunit
    !  ! Local Variables
    !  Logical :: to_close
    !  Inquire(Unit=iunit,Opened=to_close)
    !  If(to_close) Close(iunit)
    !End Subroutine dash_mpi_close_unit

  !! *******************************************
  !! ============ Private Procedures ===========
  !! *******************************************
    Function get_type(value) Result(this_datatype)
    !! Private function to retrive the type parameter for MPI procedures
    !! Only 'standard' characer should be passed to MPI procedures...
      Integer              :: this_datatype
      Class(*), Intent(In) :: value
      Select Type(value)
      Type Is(Integer)
        this_datatype = DASH_TYPE_INT
      Type Is(Real(fp64))
        this_datatype = DASH_TYPE_REAL
      Type Is(Complex(fp64))
        this_datatype = DASH_TYPE_COMP
      Type Is(Logical)
        this_datatype = DASH_TYPE_LOG
      Class Default
        Call dash_mpi_interface_stop_error(0,'get_type','invalid datatype')
      End Select
    End Function get_type

    Subroutine dash_stop
      ! Local Variables
      Integer :: max_unit = 999
      Integer :: i
      Integer :: ierr
      Logical :: to_close
      !! Try to close all the units
      Do i=1, max_unit
        Inquire(Unit=i, Opened=to_close)
        If(to_close) Close(i)
      Enddo
      Call MPI_abort(MPI_COMM_WORLD,-6,ierr)
      STOP
    End Subroutine dash_stop


















  Subroutine data_file_metadata_create(mdata,pos,size,type,data_size)
    Class(data_file_metadata), Intent(  Out) :: mdata
    Integer                  , Intent(In   ) :: pos  !! iterator value, index of k in this batch
    Integer                  , Intent(In   ) :: size !! total size of batch (ks points assigned)
    Integer                  , Intent(In   ) :: type !! this k type
    Integer                  , Intent(In   ) :: data_size !! number of basis for this k
    ! Local Variables
    Integer(fp64) :: rl         !! size of the matrix, using fp64 precision to avoid overflow
    Integer       :: rl2        !! size of one value
    Real(fp64)    :: var_real   !! just used as testing sample
    Complex(fp64) :: var_comp   !! just used as testing sample
    Select Case(type)
    Case(POINT_IS_REAL)
      Inquire(iolength=rl2) var_real
    Case(POINT_IS_COMPLEX)
      Inquire(iolength=rl2) var_comp
    Case Default
      Call dash_mpi_interface_stop_error(0,"data_file_metadata_create","incorrect point type")
    End Select
    mdata%record        = pos
    mdata%record_lenght = rl2*(data_size**2)   !! now set to the real record_lenght for this ks, shold be incresed if this is real and there are complex records
    mdata%record_size   = size
    mdata%record_type   = type
  End Subroutine data_file_metadata_create

  Subroutine data_file_open_init(file,name,mdata)
  !! Open and initialize data_file
    Class(data_file)        ,               Intent(InOut) :: file
    Character(Len=*)        ,               Intent(In   ) :: name
    Type(data_file_metadata), Dimension(:), Intent(InOut) :: mdata
    ! Local Variables
    Integer(int64) :: max_rec_len
    Integer        :: m
    If(file%initialized()) Call dash_mpi_interface_stop_error(0,"data_file_open_init","already opened")
    Allocate(file%metadata(Size(mdata)))
    max_rec_len = 0_int64
    Do m = 1, Size(mdata)
      max_rec_len = Max( mdata(m)%record_lenght, max_rec_len )
    !max_rec_len = Max( mdata(:)%record_lenght )
    Enddo
    Do m = 1, Size(mdata)
      mdata(m)%record_lenght = max_rec_len
    Enddo
    !mdata(:)%record_lenght = max_rec_len
    file%metadata = mdata
    Open( Newunit=file%unit, File='DASH-'//get_unique_id(3)//'.'//Trim(get_unique_suffix()), &
        & Action='readwrite',Form='unformatted',Access='direct',Recl=max_rec_len)
    file%status = INITIALIZED
  End Subroutine data_file_open_init

  Function data_file_initialized(file) Result(res)
  !! Check if already initialized
    Logical                      :: res
    Class(data_file), Intent(In) :: file
    res = .False.
    If(file%status.eq.INITIALIZED) res = .True.
  End Function data_file_initialized

  Subroutine data_file_close(file)
  !! Close storage file, is also possible that file is not opened by this process...
    Class(data_file), Intent(InOut) :: file
    If(file%status.eq.INITIALIZED) Then
      Close(file%unit,Status='delete')
      Deallocate(file%metadata)
      file%status = NOT_INITIALIZED
    Endif
  End Subroutine data_file_close

  Subroutine data_file_real_write_to_disk(file,buffer_data,record)
    Class(data_file),                 Intent(InOut) :: file
    Real(fp64)      , Dimension(:,:), Intent(In   ) :: buffer_data
    Integer         ,                 Intent(In   ) :: record
    integer :: rl
    If(.not.file%initialized()) Call dash_mpi_interface_stop_error(0,"data_file_real_write_to_disk","not initialized")
    Write(file%unit,Rec=record) buffer_data
  End Subroutine data_file_real_write_to_disk

  Subroutine data_file_comp_write_to_disk(file,buffer_data,record)
    Class(data_file),                 Intent(InOut) :: file
    Complex(fp64)   , Dimension(:,:), Intent(In   ) :: buffer_data
    Integer         ,                 Intent(In   ) :: record
    integer :: rl
    If(.not.file%initialized()) Call dash_mpi_interface_stop_error(0,"data_file_comp_write_to_disk","not initialized")
    Write(file%unit,Rec=record) buffer_data
  End Subroutine data_file_comp_write_to_disk

  Subroutine data_file_real_read_from_disk(file,buffer_data,record)
    Class(data_file),                 Intent(InOut) :: file
    Real(fp64)      , Dimension(:,:), Intent(  Out) :: buffer_data
    Integer         ,                 Intent(In   ) :: record
    If(.not.file%initialized()) Call dash_mpi_interface_stop_error(0,"data_file_real_read_from_disk","not initialized")
    Read(file%unit,Rec=record) buffer_data
  End Subroutine data_file_real_read_from_disk

  Subroutine data_file_comp_read_from_disk(file,buffer_data,record)
    Class(data_file),                 Intent(InOut) :: file
    Complex(fp64)   , Dimension(:,:), Intent(  Out) :: buffer_data
    Integer         ,                 Intent(In   ) :: record
    If(.not.file%initialized()) Call dash_mpi_interface_stop_error(0,"data_file_comp_read_from_disk","not initialized")
    Read(file%unit,Rec=record) buffer_data
  End Subroutine data_file_comp_read_from_disk

  Subroutine data_file_vect_real_read_from_disk(file,buffer_vect,record)
    Class(data_file),               Intent(InOut) :: file
    Real(fp64)      , Dimension(:), Intent(  Out) :: buffer_vect
    Integer         ,               Intent(In   ) :: record
    If(.not.file%initialized()) Call dash_mpi_interface_stop_error(0,"data_file_vect_real_read_from_disk","not initialized")
    Read(file%unit,Rec=record) buffer_vect
  End Subroutine data_file_vect_real_read_from_disk

  Subroutine data_file_transfer_ownership(file,file2)
    Class(data_file), Intent(InOut) :: file
    Type(data_file) , Intent(InOut) :: file2
    Call file%close()
    !! From now on only relevant processes should trasfer the file....
    If(file2%status.eq.INITIALIZED) Then
      file%unit     = file2%unit
      file%metadata = file2%metadata
      file%status = INITIALIZED
      ! ---------------------
      file2%unit  = INVALID
      Deallocate(file2%metadata)
      file2%status = NOT_INITIALIZED
    Endif
  End Subroutine data_file_transfer_ownership

  Function get_unique_suffix() Result(suffix)
    Character(Len=10) :: suffix
    write(suffix,"(A3,I0)") 'rnk', world%rank
  End Function get_unique_suffix

  Function get_unique_id(len) Result(id)
    Integer           , Intent(In   ) :: len
    Character(Len=len)                :: id
    ! Local Variables
    Character(Len=36), Parameter :: charset = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Integer          , Save      :: current_index = 1  !! keeps track of current count
    Integer                      :: temp_index
    Integer                      :: i, digit
    Character(Len=len)           :: res
    res = Repeat('0',len)
    temp_index = current_index
    Do i = len, 1, -1
      digit = Mod(temp_index,36)
      res(i:i) = charset(digit+1:digit+1)
      temp_index = temp_index / 36
      If(temp_index.eq.0) Exit
    Enddo
    id = res
    current_index = current_index + 1
  End Function get_unique_id














  !! **********************************************
  !! ============ Type-bound Procedures ===========
  !! **********************************************
    Subroutine proc_info_create_from_comm(R,old_comm)
    !! Generate the generic proc_info object and populate its components
      Class(proc_info), Intent(InOut) :: R
      Integer         , Intent(In   ) :: old_comm
      ! Local Variables
      Integer :: ierr
      Integer :: my_rank
      Integer :: my_size
      !! First check if already initialized
      If(R%status.eq.INITIALIZED)Then
        Call dash_mpi_interface_stop_error(0,'proc_info_create','already initialized')
      Endif
      !! Setting new communicator
      Call MPI_comm_dup(old_comm,R%comm,ierr)
      !! Retrive standard parallelization info
      Call MPI_comm_rank(R%comm,my_rank,ierr)
      Call MPI_comm_size(R%comm,my_size,ierr)
      R%rank = my_rank
      R%size = my_size
      !! Update status
      R%status = INITIALIZED
    End Subroutine proc_info_create_from_comm

    Function proc_info_is_initialized(R) Result(is_initialized)
      Logical                      :: is_initialized
      Class(proc_info), Intent(In) :: R
      Select Case(R%status)
      Case(INITIALIZED)
        is_initialized = .True.
      Case(NOT_INITIALIZED)
        is_initialized = .False.
        Call dash_mpi_interface_stop_error(1,'proc_info_is_initialized','not initialized')
      Case Default
        is_initialized = .False.
        Call dash_mpi_interface_stop_error(0,'proc_info_is_initialized','invalid status')
      End Select
    End Function proc_info_is_initialized

    ! --- pointers for paralleliztions like in PCRYSTAL ---
    Subroutine proc_info_mp_reset(R)
    !! Set initial value of my_point
    !! Inverse numbering as in PCRYSTAL
    !! From 1 to size
      Class(proc_info), Intent(InOut) :: R
      If(R%size.le.0) Call dash_mpi_interface_stop_error(0,'proc_info_mp_reset','communicator size .le. 0')
      R%my_point = R%size - R%rank
    End Subroutine proc_info_mp_reset

    Subroutine proc_info_mp_next(R)
    !! Increment my_point value
      Class(proc_info), Intent(InOut) :: R
      If(R%size.le.0) Return
      R%my_point = R%my_point + R%size
    End Subroutine proc_info_mp_next

    ! --- MPI operations ---
    Subroutine proc_info_sync(R)
    !! Sync all ranks which belong to this communicator
      Class(proc_info), Intent(In) :: R
      ! Local Variables
      Integer :: ierr
      If(R%comm.ne.MPI_COMM_NULL) Call MPI_barrier(R%comm,ierr)
    End Subroutine proc_info_sync

    Subroutine proc_info_max(R,value)
    !! Compute maximum value between all processes in this communicator
    !! Overwrite value
      Class(proc_info), Intent(In   ) :: R
      Class(*)        , Intent(InOut) :: value
      ! Local Variables
      Integer :: ierr
      If(R%comm.ne.MPI_COMM_NULL) Then
        Call MPI_allreduce(MPI_IN_PLACE,value,1,get_type(value),MPI_MAX,R%comm,ierr)
      Endif
    End Subroutine proc_info_max

    Subroutine proc_info_sum(R,value)
    !! Compute summation value between all processes in this communicator
    !! Overwrite value
      Class(proc_info), Intent(In   ) :: R
      Class(*)        , Intent(InOut) :: value
      ! Local Variables
      Integer :: ierr
      If(R%comm.ne.MPI_COMM_NULL) Then
        Call MPI_allreduce(MPI_IN_PLACE,value,1,get_type(value),MPI_SUM,R%comm,ierr)
      Endif
    End Subroutine proc_info_sum
    Subroutine proc_info_sum_r1(R,value)
    !! Compute summation value between all processes in this communicator
    !! Overwrite value
    !! value is a rank 1, the polimorfism using select rank is not feasible due
    !! to intel fortran compiler internal bug...
      Class(proc_info),               Intent(In   ) :: R
      Class(*)        , Dimension(:), Intent(InOut) :: value
      ! Local Variables
      Integer :: ierr
      If(R%comm.ne.MPI_COMM_NULL) Then
        Call MPI_allreduce(MPI_IN_PLACE,value,Size(value),get_type(value(1)),MPI_SUM,R%comm,ierr)
      Endif
    End Subroutine proc_info_sum_r1


    Subroutine proc_info_broadcast_r1(R,value,root)
    !! ============================================================
    !! Broadcast 1D vector value from root process to all processes 
    !! in R%comm (assuming root is a valid rank for R%comm)
    !!  * Every process should have value allocated with same size
    !! ============================================================
    !! WARNING dont know if polymorphism is well digested here...
      Class(proc_info),               Intent(In   ) :: R
      Class(*)        , Dimension(:), Intent(InOut) :: value
      Integer         ,               Intent(In   ) :: root
      ! Local Variables
      Integer :: ierr
      If(R%comm.ne.MPI_COMM_NULL) Then
        CALL MPI_BCAST(value,Size(value),get_type(value(1)),root,R%comm,ierr)
      Endif
    End Subroutine proc_info_broadcast_r1

    Subroutine proc_info_gather(R,send_value,receve_value)
    !! Get all send_values to rank 0 process, automatic allocation of receve_value
      Class(proc_info),                            Intent(In   ) :: R
      Integer         ,                            Intent(In   ) :: send_value
      Integer         , Dimension(:), Allocatable, Intent(InOut) :: receve_value
      ! Local Variables
      Integer :: ierr
      Integer :: n_conunt
      n_conunt = R%size
      If(Allocated(receve_value)) Deallocate(receve_value)
      If(R%rank.eq.0) Allocate(receve_value(n_conunt))
      If(R%comm.ne.MPI_COMM_NULL) Then
        Call MPI_gather( send_value,1,MPI_INTEGER,receve_value,1, &
                       & MPI_INTEGER,0,R%comm,ierr )
      Endif
    End Subroutine proc_info_gather
    Subroutine proc_info_gather_r1(R,send_value,receve_value)
    !! Get all send_values to rank 0 process, automatic allocation of receve_value
    !! For 1D array to send, should all be of same size
      Class(proc_info),                            Intent(In   ) :: R
      Integer         , Dimension(:),              Intent(In   ) :: send_value
      Integer         , Dimension(:), Allocatable, Intent(InOut) :: receve_value
      ! Local Variables
      Integer :: ierr
      Integer :: n_conunt
      Integer :: send_size, max_send_size
      n_conunt = R%size
      send_size = Size(send_value, Dim=1)
      max_send_size = send_size
      Call R%max(max_send_size)
      If(send_size.ne.max_send_size) Call dash_mpi_interface_stop_error(0,"proc_info_gather_r1","send buffers not same size")
      If(Allocated(receve_value)) Deallocate(receve_value)
      If(R%rank.eq.0) Allocate(receve_value( n_conunt * send_size ))
      If(R%comm.ne.MPI_COMM_NULL) Then
        Call MPI_gather( send_value,send_size,MPI_INTEGER,receve_value,send_size, &
                       & MPI_INTEGER,0,R%comm,ierr )
      Endif
    End Subroutine proc_info_gather_r1
    Subroutine proc_info_gather_real(R,send_value,receve_value)
    !! Get all send_values to rank 0 process, automatic allocation of receve_value
      Class(proc_info),                            Intent(In   ) :: R
      Real(fp64)      ,                            Intent(In   ) :: send_value
      Real(fp64)      , Dimension(:), Allocatable, Intent(InOut) :: receve_value
      ! Local Variables
      Integer :: ierr
      Integer :: n_conunt
      n_conunt = R%size
      If(Allocated(receve_value)) Deallocate(receve_value)
      If(R%rank.eq.0) Allocate(receve_value(n_conunt))
      If(R%comm.ne.MPI_COMM_NULL) Then
        Call MPI_gather( send_value,1,MPI_DOUBLE_PRECISION,receve_value,1, &
                       & MPI_DOUBLE_PRECISION,0,R%comm,ierr )
      Endif
    End Subroutine proc_info_gather_real

    ! --- extract info ---
    Subroutine proc_info_get(R,comm,rank,size,status,my_point)
    !! Extract informations about this proc_info object
      Class(proc_info), Intent(In   )           :: R
      Integer         , Intent(  Out), Optional :: comm
      Integer         , Intent(  Out), Optional :: rank
      Integer         , Intent(  Out), Optional :: size
      Integer         , Intent(  Out), Optional :: status
      Integer         , Intent(  Out), Optional :: my_point
      If(Present(comm    )) comm     = R%comm
      If(Present(rank    )) rank     = R%rank
      If(Present(size    )) size     = R%size
      If(Present(status  )) status   = R%status
      If(Present(my_point)) my_point = R%my_point
    End Subroutine proc_info_get

End Module interface_MPI
