Module module_arrays
!===========================================================
!    @author: Giacomo Ambrogio       date: Feb 2024
!
!    DASH - Device Accelerated Solutions for HPC
!
!    K point array main class
!
!    This module contains the outermost data structures
!    for dealing with k points
!
!                                    date: May 2025
!    Now using pointers for dealing with ks points info
!===========================================================
Use base_numbers
Use module_points
Use interface_MPI, Only: data_file, data_file_metadata

  Implicit None
  Private

  Integer, Parameter, Public :: N_IRREPS_NOT_DEFINED = -1
  ! Should implement a get_irreps to get the irreps from a ks_array
  ! based on the iterator value, shoyuld call the get_irreps from ks_point

  Integer, Public, Protected :: memory_execution_mode
  Integer, Public, Protected :: memory_storage_mode

  Type, Public :: info_ks_point
    Integer   ,                            Private :: k_index
    Integer   ,                            Private :: s_index
    Integer   ,                            Private :: ks_type
    Integer   ,                            Private :: ks_size
    Integer   ,                            Private :: irrep      = INVALID    !! Not yet implemented
    Integer   , Dimension(3),              Private :: coord
    Real(fp64),                            Private :: weight     = INVALID    !! Weight of the k point for the integration over k. sum( 1/nk ) * mvf * comp|real
    Integer   ,                            Private :: occ_states
    Real(fp64), Dimension(:), Allocatable, Private :: occ_bands
    Contains
    Procedure, Public :: create    => info_ks_point_create
    Procedure, Public :: update    => info_ks_point_update
    Procedure, Public :: exist     => info_ks_point_exist
    Procedure, Public :: get_type  => info_ks_point_get_type
    Procedure, Public :: get_coord => info_ks_point_get_coord
    Procedure, Public :: get_spin  => info_ks_point_get_spin
    Procedure, Public :: get_k     => info_ks_point_get_k
    Procedure, Public :: get_occ   => info_ks_point_get_occ
    Procedure, Public :: get_size  => info_ks_point_get_size
    Procedure, Public :: is_gamma  => info_ks_point_is_gamma
  End Type


  Type, Public :: ks_point_data
    Type(ks_point)     ,          Private :: data
    Type(info_ks_point), Pointer, Private :: info   => Null()
    Integer            ,          Private :: status =  POINT_STATUS_DEFAULT
  End Type

  Type, Public :: ks_point_data_1D
    Type(ks_point_1D)  ,          Private :: data
    Type(info_ks_point), Pointer, Private :: info   => Null()
    Integer            ,          Private :: status =  POINT_STATUS_DEFAULT
  End Type


  Type, Public :: ks_array_1D
  !! Type to contains 1D arrays, follows the same structure as the main array class
    Type(ks_point_data_1D), Dimension(:), Allocatable, Private :: data
    Type(info_ks_point)   , Dimension(:), Allocatable, Private :: info
    Integer               ,                            Private :: status          =  POINT_STATUS_DEFAULT
    Type(info_ks_point)   , Dimension(:), Pointer    , Private :: all_points_info => Null()
    Integer               ,                            Private :: iterator_value  =  INVALID
    Contains
    Procedure, Public :: create         => ks_array_1D_create
    Procedure, Public :: destroy        => ks_array_1D_destroy
    Procedure, Public :: get_raw        => ks_array_1D_get_raw
    Procedure, Public :: set_raw        => ks_array_1D_set_raw
    Procedure, Public :: shift          => ks_array_1D_shift
    Procedure, Public :: to_ene         => ks_array_1D_to_ene
    Procedure, Public :: print          => ks_array_1D_print
    Procedure, Public :: compute_fermi  => ks_array_1D_compute_fermi
    Procedure, Public :: sort           => ks_array_1D_sort
    !Procedure, Public :: from_ene =>  !! TO DO
    Procedure, Public :: iterator_init  => ks_array_1D_iterator_init
    Procedure, Public :: iterator_reset => ks_array_1D_iterator_reset
    Procedure, Public :: iterator_next  => ks_array_1D_iterator_next
  End Type

  Type, Public :: ks_array
    Type(ks_point_data), Dimension(:), Allocatable, Private :: ks_points
    Integer                                       , Private :: n_ks_points      =  0
    Integer            , Dimension(2)             , Private :: n_ks_points_type =  (0,0)
    Type(info_ks_point), Dimension(:), Pointer    , Private :: all_points_info  => Null()
    Integer            , Dimension(:), Allocatable, Private :: my_indeces       !! This is the list of owned ks. If no ks assigned then %my_indeces is not allocated
    Integer                                       , Private :: iterator_value   = INVALID
    Type(data_file)                               , Private :: storage
    Integer                                       , Private :: status           = STATUS_INVALID
    Contains
      !! --- LA Operations ---
    Procedure, Public :: initialize_ortog_mat => ks_array_initialize_ortogonalization_matrices
    Procedure, Public :: multiply             => ks_array_multiply
    Procedure, Public :: similarity           => ks_array_similarity
    Procedure, Public :: diag                 => ks_array_diag
    Procedure, Public :: copy                 => ks_array_copy
      !! --- Iterative Operations ---
    Procedure, Public :: shift_diagonal       => ks_array_shift_diagonal
      !! --- Utilities ---
    Generic  , Public :: create               => ks_array_create
    Generic  , Public :: create               => ks_array_create_vary
    Procedure, Public :: destroy              => ks_array_destroy
    Procedure, Public :: global_store         => global_store_ks_array
    Procedure, Public :: global_free          => global_free_ks_array
    Procedure, Public :: to_store             => ks_array_to_store
    Procedure, Public :: from_store           => ks_array_from_store
    Procedure, Public :: allocated            => ks_array_allocated
    Procedure, Public :: initialized          => ks_array_initialized
    Procedure, Public :: size                 => ks_array_get_number_of_ks_points
    Procedure, Public :: sizeof               => ks_array_get_device_mem_estimation
    Procedure, Public :: print_mem            => ks_array_print_memory_report
    Procedure, Public :: set_indexing         => ks_array_set_indexing_from_array
    Procedure, Public :: print_dist           => ks_array_print_distribution
    Generic  , Public :: set_raw              => set_raw_real, set_raw_comp
    Generic  , Public :: get_raw              => get_raw_real, get_raw_comp, get_raw_vect
    Procedure, Public :: compute_tester       => ks_array_compute_tester
    Procedure, Public :: iterator_init        => ks_array_iterator_init
    Procedure, Public :: iterator_next        => ks_array_iterator_next
    Procedure, Public :: iterator_reset       => ks_array_iterator_reset
    Procedure, Public :: compute_density_mat  => ks_array_compute_density
    ! Private Implementations
    Procedure, Private :: ks_array_create
    Procedure, Private :: ks_array_create_vary
    Procedure, Private :: set_raw_real                    => ks_array_set_raw_real
    Procedure, Private :: set_raw_comp                    => ks_array_set_raw_comp
    Procedure, Private :: get_raw_real                    => ks_array_get_raw_real
    Procedure, Private :: get_raw_vect                    => ks_array_get_raw_vector
    Procedure, Private :: get_raw_comp                    => ks_array_get_raw_comp
    Procedure, Private :: load_data                       => ks_array_load_data
    Procedure, Private :: dump_data                       => ks_array_dump_data
  End Type

  Public :: dash_arrays_set_exec_mem_option
  Public :: dash_arrays_set_storage_mem_option

  Public :: get_maximum_k_type

  !! ========================= STORAGE =================================
  Type(ks_array)     ,                                    Public :: base
  Type(ks_array)     ,                                    Public :: sk_array
  Type(ks_array)     ,                                    Public :: hk_array
  Type(ks_array)     ,                                    Public :: fk_array
  Type(ks_array_1D)  ,                                    Public :: evals_array
  Type(info_ks_point), Dimension(:), Allocatable, Target, Public :: all_ks_points_info_target
  !! ========================= STORAGE =================================

  Contains

    Subroutine dash_arrays_set_exec_mem_option(mode)
      Integer, Intent(In) :: mode !! already validated in dash_module
      memory_execution_mode = mode
    End Subroutine dash_arrays_set_exec_mem_option

    Subroutine dash_arrays_set_storage_mem_option(mode)
      Integer, Intent(In) :: mode !! already validated in dash_module
      memory_storage_mode = mode
    End Subroutine dash_arrays_set_storage_mem_option

!!#############################################################################
!!        info_ks_point type-bound procedures
!!#############################################################################
    Subroutine info_ks_point_create(I,k_index,s_index,ks_type,ks_size,n_of_irreps,n_occ_ao,k_coords,weight)
    !! Reset all info of this ks, and sets new if present
    !! Used to generete the first set of info_ks_point
    Use dash_utils, Only: dash_error
      Class(info_ks_point),               Intent(InOut)           :: I
      Integer             ,               Intent(In   ), Optional :: k_index
      Integer             ,               Intent(In   ), Optional :: s_index
      Integer             ,               Intent(In   ), Optional :: ks_type
      Integer             ,               Intent(In   ), Optional :: ks_size
      Integer             ,               Intent(In   ), Optional :: n_of_irreps
      Integer             ,               Intent(In   ), Optional :: n_occ_ao
      Integer             , Dimension(3), Intent(In   ), Optional :: k_coords
      Real(fp64)          ,               Intent(In   ), Optional :: weight
      If(Present(n_of_irreps).and.(n_of_irreps.gt.1)) Then
        Call dash_error(0,"info_ks_point_create","number of irreps > 1")
      Endif
      I%k_index    = POINT_NOT_EXIST
      I%s_index    = SPIN_ALPHA
      I%ks_type    = POINT_NOT_EXIST
      I%ks_size    = 1
      I%irrep      = 1
      I%coord      = [0,0,0]
      I%weight     = INVALID
      I%occ_states = 0
      If(Present(k_index    )) I%k_index    = k_index
      If(Present(s_index    )) I%s_index    = s_index
      If(Present(ks_type    )) I%ks_type    = ks_type
      If(Present(ks_size    )) I%ks_size    = ks_size
      If(Present(n_of_irreps)) I%irrep      = n_of_irreps
      If(Present(k_coords   )) I%coord      = k_coords
      If(Present(weight     )) I%weight     = weight
      If(Present(n_occ_ao   )) I%occ_states = n_occ_ao
      If(Allocated(I%occ_bands)) Deallocate(I%occ_bands)
      Allocate(I%occ_bands(I%ks_size))
      I%occ_bands = -9.0_fp64
    End Subroutine info_ks_point_create

    Subroutine info_ks_point_update(I,ks_size,nocc,occupations)
    Use dash_utils, Only: dash_error
      Class(info_ks_point) ,               Intent(InOut)           :: I
      Integer              ,               Intent(In   ), Optional :: ks_size
      Integer              ,               Intent(In   ), Optional :: nocc
      Real(fp64)           , Dimension(:), Intent(In   ), Optional :: occupations
      ! Local Valiable
      Character(Len=12) :: str
      !! --- Update Size ---
      If(Present(ks_size)) Then
        I%ks_size = ks_size
        If(Allocated(I%occ_bands)) Deallocate(I%occ_bands)
        Allocate(I%occ_bands(ks_size))
        I%occ_bands = -9.0_fp64
        Write(str,"(3(I4))") I%coord(1),I%coord(2),I%coord(3)
        Call dash_error(1,"info_ks_point_update","Updating size of k point "//str)
      Endif
      !! --- Update Occupation ---
      If(Present(nocc).and.Present(occupations)) Then
        I%occ_states = nocc
        If(Size(occupations).ne.I%ks_size) Call dash_error(0,"info_ks_point_update","wrong occupations size")
        I%occ_bands = occupations
        Write(str,"(2(I6))") nocc, Count(Abs(occupations).gt.1.0E-11_fp64)
        If( nocc.ne.Count(Abs(occupations).gt.1.0E-11_fp64) ) Call dash_error(1,"info_ks_point_update","nocc and occupations not coherent "//str)
      Else
        Call dash_error(0,"info_ks_point_update","need both occupation number and bands")
      Endif

    !    !! -----------> NO FRACTIONAL OCCUPATION (ONLY FOR MOLECULES)
    !    !If(lnofoccu) Then
    !    !  alfa = 0._fp64
    !    !  Do ijq = 1, nband
    !    !    If(jqoccs(ijq).ne.0) Then
    !    !      alfa(ijq) = 1._fp64
    !    !    Endif
    !    !  Enddo
    !    !Endif
    !    !! <-----------

    End Subroutine info_ks_point_update

    Logical Function info_ks_point_exist(I)
    Use dash_utils, Only: dash_error
      Class(info_ks_point), Intent(InOut) :: I
      Selectcase(I%ks_type)
      Case(POINT_NOT_EXIST)
        info_ks_point_exist = .False.
      Case(POINT_IS_REAL,POINT_IS_COMPLEX)
        info_ks_point_exist = .True.
      Case Default
        info_ks_point_exist = .False.
        Call dash_error(0,'info_ks_point_exist','invalid k type')
      Endselect
    End Function info_ks_point_exist

    Integer Function info_ks_point_get_type(I)
      Class(info_ks_point), Intent(In) :: I
      info_ks_point_get_type = I%ks_type
    End Function info_ks_point_get_type

    Integer Function info_ks_point_get_coord(I,d)
      Class(info_ks_point), Intent(In) :: I
      Integer             , Intent(In) :: d
      info_ks_point_get_coord = I%coord(d)
    End Function info_ks_point_get_coord

    Integer Function info_ks_point_get_spin(I)
      Class(info_ks_point), Intent(In) :: I
      info_ks_point_get_spin = I%s_index
    End Function info_ks_point_get_spin

    Integer Function info_ks_point_get_k(I)
      Class(info_ks_point), Intent(In) :: I
      info_ks_point_get_k = I%k_index
    End Function info_ks_point_get_k

    Integer Function info_ks_point_get_occ(I)
      Class(info_ks_point), Intent(In) :: I
      info_ks_point_get_occ = I%occ_states
    End Function info_ks_point_get_occ

    Integer Function info_ks_point_get_size(I)
      Class(info_ks_point), Intent(In) :: I
      info_ks_point_get_size = I%ks_size
    End Function info_ks_point_get_size

    Logical Function info_ks_point_is_gamma(I)
      Class(info_ks_point), Intent(In) :: I
      info_ks_point_is_gamma = I%coord(1).eq.0 .and. &
                             & I%coord(2).eq.0 .and. &
                             & I%coord(3).eq.0
    End Function info_ks_point_is_gamma

    Integer Function get_maximum_k_type(I_array)
      Type(info_ks_point), Dimension(:), Intent(In) :: I_array
      ! Local Variables
      Integer :: i
      Logical :: comp_present
      comp_present = .False.
      Do i = 1, Size(I_array)
        If(I_array(i)%ks_type.eq.POINT_IS_COMPLEX) comp_present = .True.
      Enddo
      get_maximum_k_type = Merge(POINT_IS_COMPLEX,POINT_IS_REAL,comp_present)
    End Function get_maximum_k_type

!!#############################################################################



!!#############################################################################
!!        ks_array type-bound procedures
!!#############################################################################
!! ========== Public Implementations ==========
!! --- Utilities ---
    Function ks_array_allocated(A) Result(res)
    !! Returns the number of allocated (host/device) points
      Integer                        :: res
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Type(info_ks_point) :: this_ks
      res = 0
      If(.not.A%initialized()) Return
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        If(A%ks_points(A%iterator_value)%data%allocated()) res = res + 1
        this_ks = A%iterator_next()
      Enddo
    End Function ks_array_allocated

    Function ks_array_initialized(A) Result(res)
    !! Returns logical, ks point vector is allocated?
    !! If array has no ks then returns always False
      Logical :: res
      Class(ks_array), Intent(In) :: A
      res = Allocated(A%ks_points)
    End Function ks_array_initialized

    Function ks_array_get_number_of_ks_points(A) Result(res)
    Use dash_utils, Only: dash_error
      Integer                     :: res
      Class(ks_array), Intent(In) :: A
      res = A%n_ks_points
      !! just a check
      If(A%initialized()) Then
        If(Size(A%ks_points).ne.res) Call dash_error(0,'ks_array_get_number_of_ks_points','size of ks_array not coherent')
      Else
        If(res.ne.0) Call dash_error(0,'ks_array_get_number_of_ks_points','size of ks_array not coherent')
      Endif
    End Function ks_array_get_number_of_ks_points

    Function ks_array_get_device_mem_estimation(A) Result(mem)
    !! ===========================================================
    !! Returns the estimated memory (MB) required to store all the
    !! points assigned to this array in device memory
    !!
    !! --- Do not take into account diagonalization workspace ---
    !! ---                  and ft_map memory                 ---
    !!
    !! The formula used is:
    !!    mem = ( 3 * \sum_{ks} ( byte_size_{R/C} * AOs^2 ) + 
    !!                \sum_{ks} ( byte_size_{R} * AOs )   ) /
    !!                          1024^2
    !!
    !!    Where byte_size_{R/C} is the bytes used to store 
    !!    real or complex values
    !!    The 3 take int account that for mat mult we need
    !!    3 arrays, and the second line is the one for the
    !!    eigenvalues
    !! ===========================================================
    Use dash_utils
      Real(fp64)                     :: mem
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Integer               :: i
      Type(info_ks_point)   :: this_ks
      Real(fp64)            :: real_sample
      Complex(fp64)         :: complex_sample
      Real(fp64), Parameter :: R_SIZE = Real(Storage_size(real_sample),kind=fp64) / 8.0_fp64
      Real(fp64), Parameter :: C_SIZE = Real(Storage_size(complex_sample),kind=fp64) / 8.0_fp64
      mem = 0.0_fp64
      Do i = 1, A%n_ks_points
        If(.not.Allocated(A%my_indeces)) Call dash_error(0,"ks_array_get_device_mem_estimation","indeces not allocated")
        this_ks = A%all_points_info(A%my_indeces(i))
        Select Case(this_ks%ks_type)
        Case(POINT_IS_REAL)
          mem = mem +  (   3 * ( R_SIZE * (this_ks%ks_size**2) ) / (1024._fp64**2)   )
        Case(POINT_IS_COMPLEX)
          mem = mem +  (   3 * ( C_SIZE * (this_ks%ks_size**2) ) / (1024._fp64**2)   )
        End Select
        mem = mem +  (   ( R_SIZE * this_ks%ks_size ) / (1024._fp64**2)   )
      Enddo
    End Function ks_array_get_device_mem_estimation

    Subroutine ks_array_print_memory_report(A,mem_lim,unit)
    !! ===========================================================
    !! Print memory report for different processes
    !!
    !! This procedure is not very clean... is quite dependent
    !! On the way we handle memory limit in dash_module, so need
    !! to be extra carefull
    !! ===========================================================
    Use dash_utils
    Use interface_MPI
      Class(ks_array), Intent(InOut) :: A
      Real(fp64)     , Intent(In   ) :: mem_lim
      Integer        , Intent(In   ) :: unit
      ! Local Variables
      Integer                               :: proc,max_procs
      Real(fp64)                            :: est_mem
      Real(fp64), Dimension(:), Allocatable :: buffer_mem
      Write(unit,"(1X,79('*'))")
      Write(unit,"(24X,A)") '*** GCRYSTAL MEMORY USAGE REPORT ***'

      Write(unit,"(1X,A)") "KS ARRAY STORAGE LOCATION:"
      Select Case(memory_storage_mode)
      Case(MEMORY_STORAGE_DEVICE)
        Write(unit,"(3X,'-',1X,A)") "DURING NON-LA STEPS: DEVICE"
      Case(MEMORY_STORAGE_HOST)
        Write(unit,"(3X,'-',1X,A)") "DURING NON-LA STEPS: HOST"
      Case(MEMORY_STORAGE_DISK)
        Write(unit,"(3X,'-',1X,A)") "DURING NON-LA STEPS: DISK"
      End Select
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DEVICE)
        Write(unit,"(3X,'-',1X,A)") "DURING LA STEPS    : DEVICE"
        Write(unit,"(1X)")
        If(mem_lim.gt.0) Then
          Write(unit,"(1X,A,F15.3,A)") "MEMORY LIMIT:",mem_lim," MB"
        Else
          Write(unit,"(1X,A)") "NO MEMORY LIMIT SET"
        Endif
      Case(MEMORY_EXECUTION_HOST)
        Write(unit,"(3X,'-',1X,A)") "DURING LA STEPS    : HOST"
      Case(MEMORY_EXECUTION_DISK)
        Write(unit,"(3X,'-',1X,A)") "DURING LA STEPS    : DISK"
      End Select
      Write(unit,"(1X)")

      Write(unit,"(1X,A)") "PER-PROCESS MEMORY USAGE:"
      Call world%get(size=max_procs)
      est_mem = A%sizeof()
      Call world%gather(est_mem,buffer_mem)
      If(Allocated(buffer_mem)) Then
        Do proc = 1, max_procs
          If(buffer_mem(proc).le.0) Cycle
          Write(unit,"(1X,'P',I5,': ')", Advance='no') proc - 1
          Write(unit,"(F15.2,A)", Advance='no') buffer_mem(proc),' MB'
          If( mem_lim.gt.0 .and. &
            & buffer_mem(proc).ge.mem_lim .and. &
            & memory_execution_mode.eq.MEMORY_EXECUTION_DEVICE ) Then
            Select Case(memory_storage_mode)
            Case(MEMORY_STORAGE_DISK,MEMORY_STORAGE_DEVICE)
              Write(unit,"(1X,A)") "  ===> EXCEEDS LIMIT — STORING ON DISK"
            Case(MEMORY_STORAGE_HOST)
              Write(unit,"(1X,A)") "  ===> EXCEEDS LIMIT — STORING ON HOST"
            End Select
          Else
            Write(unit,*)
          Endif
        Enddo
        Deallocate(buffer_mem)
      Endif
      Write(unit,"(1X,79('*'))")
    End Subroutine ks_array_print_memory_report

    Subroutine global_store_ks_array(A,xname)
    !! ============================================================
    !! Procedure to load the array in memory, different behaviours
    !! based on the memory exacution mode selected:
    !! --> DEVICE: allocate all the matrices in memory
    !! --> HOST  : allocate all the matrices in host memory 
    !! --> DISK  : * create memtadata array (description of the array
    !!               of ks points that will be stored on disk).
    !!             This array is just used to set some information
    !!             for the actual initialization of storage file.
    !!             * initialization of storage file:
    !!               - Setting max record lenght from all metadata
    !!               - Open new unit, file_unique_id 
    !!               - file%status = INITIALIZED
    !!
    !! Array needs to be initialized ( Call A%crate )
    !!
    !! Set array status accordingly
    !!
    !! At the moment we don't have any data actually stored, 
    !! but memory is reserved
    !! ============================================================
    Use dash_utils
      Class(ks_array) , Intent(InOut)           :: A
      Character(Len=*), Intent(In)   , Optional :: xname
      ! Local variables
      Integer                                             :: ks
      Integer                 , Dimension(:), Allocatable :: types
      Type(info_ks_point)                                 :: this_ks
      Character(Len=7)                                    :: name
      Type(data_file_metadata), Dimension(:), Allocatable :: all_metatdata
      If(A%size().eq.0) Then
        A%status = STATUS_INIT
        Return
      Endif
      If(.not.A%initialized()) Call dash_error(0,"global_store_ks_array","array not initialized")
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DEVICE)
        this_ks = A%iterator_init()
        Do While (this_ks%exist())
          If(.not.A%ks_points(A%iterator_value)%data%allocated()) Then
            Call A%ks_points(A%iterator_value)%data%store(1,this_ks%ks_size)
          Endif
          this_ks = A%iterator_next()
        Enddo
        A%status = STATUS_DEVICE
      Case(MEMORY_EXECUTION_DISK)
        Allocate(all_metatdata(A%n_ks_points))
        this_ks = A%iterator_init()
        Do While (this_ks%exist())
          Call all_metatdata(A%iterator_value)%create( A%iterator_value,A%n_ks_points,this_ks%ks_type,this_ks%ks_size )
          this_ks = A%iterator_next()
        Enddo
        name = 'unknown'
        If(Present(xname)) name = xname
        Call A%storage%init(name,all_metatdata)
        Deallocate(all_metatdata)
        A%status = STATUS_DISK
      Case(MEMORY_EXECUTION_HOST)
        this_ks = A%iterator_init()
        Do While (this_ks%exist())
          Call A%ks_points(A%iterator_value)%data%host_store(1,this_ks%ks_size)
          this_ks = A%iterator_next()
        Enddo
        A%status = STATUS_HOST
      End Select
    End Subroutine global_store_ks_array

    Subroutine global_free_ks_array(A)
    !! ============================================================
    !! Procedure to free memory space of the array, different 
    !! behaviours based on the memory storage mode selected:
    !! --> DEVICE: free all the matrices from device memory
    !! --> HOST  : free all the matrices from host memory
    !! --> DISK  : close unit storage file, and delete it
    !!             (we don’t want to garbage-fill the working dir)
    !!
    !! Remove status
    !!
    !! ============================================================
      Class(ks_array), Intent(InOut) :: A
      ! Local variables
      Type(info_ks_point) :: this_ks
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DEVICE)
        this_ks = A%iterator_init()
        Do While (this_ks%exist())
          If(A%ks_points(A%iterator_value)%data%allocated()) Then
            Call A%ks_points(A%iterator_value)%data%free()
          Endif
          this_ks = A%iterator_next()
        Enddo
      Case(MEMORY_EXECUTION_DISK)
        Call A%storage%close()
      Case(MEMORY_EXECUTION_HOST)
        this_ks = A%iterator_init()
        Do While (this_ks%exist())
          Call A%ks_points(A%iterator_value)%data%host_free()
          this_ks = A%iterator_next()
        Enddo
      End Select
      A%status = STATUS_INIT
    End Subroutine global_free_ks_array

    Subroutine ks_array_to_store(A)
    !! ============================================================
    !! Procedure used to store "away" the ks_array
    !! Useful when need to free space on device/host for other
    !! calculations (maybe integrals?)
    !! 
    !! Usage: after the array has been cretated and %global_store-d
    !!        Call A%to_store() will move the array data to storage
    !!        location
    !!        Before using the array, call A%from_store()
    !!
    !! This procedure is general for the whole array, to handle 
    !! single ks points use %load_data/%dump_data
    !!
    !! This is ment to be used before (and %from_store after) the
    !! LA section in the SCF (fdik + fermi + pdig)
    !!
    !! This procedure is not super-clean, but should work fine
    !! Possible improvment on memory movment when buffers are used
    !! ============================================================
    Use dash_utils
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Type(info_ks_point)                                   :: this_ks
      Type(data_file_metadata), Dimension(:)  , Allocatable :: all_metatdata
      Real(fp64)              , Dimension(:,:), Allocatable :: real_buffer
      Complex(fp64)           , Dimension(:,:), Allocatable :: comp_buffer
      If(A%size().eq.0) Return
      If(.not.A%initialized()) Call dash_error(0,"ks_array_to_store","array not initialized")
      !! If execution is equal to store nothing to do 
      If(memory_execution_mode.eq.memory_storage_mode) Return

      Select Case(memory_execution_mode)
      !! **********************************************************
      Case(MEMORY_EXECUTION_DEVICE)
      !! **********************************************************
        If(A%status.ne.STATUS_DEVICE) Call dash_error(0,"ks_array_to_store","data not on device")
        Selectcase(memory_storage_mode)
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_HOST)
          this_ks = A%iterator_init()
          Do While(this_ks%exist())
            Call A%ks_points(A%iterator_value)%data%host_store(1,this_ks%ks_size)
            Call A%ks_points(A%iterator_value)%data%dev_to_host()
            Call A%ks_points(A%iterator_value)%data%free()
            this_ks = A%iterator_next()
          Enddo
          A%status = STATUS_HOST
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_DISK)
          If(.not.A%storage%initialized()) Then
            Allocate(all_metatdata(A%n_ks_points))
            this_ks = A%iterator_init()
            Do While (this_ks%exist())
              Call all_metatdata(A%iterator_value)%create( A%iterator_value,A%n_ks_points,this_ks%ks_type,this_ks%ks_size )
              this_ks = A%iterator_next()
            Enddo
            Call A%storage%init('tmp',all_metatdata)
            Deallocate(all_metatdata)
          Endif
          this_ks = A%iterator_init()
          Do While (this_ks%exist())
            Select Case(this_ks%ks_type)
            Case(POINT_IS_REAL)
              Call A%ks_points(A%iterator_value)%data%get_raw(real_buffer)
              Call A%storage%data_to_disk(real_buffer,A%iterator_value)
              Deallocate(real_buffer)
            Case(POINT_IS_COMPLEX)
              Call A%ks_points(A%iterator_value)%data%get_raw(comp_buffer)
              Call A%storage%data_to_disk(comp_buffer,A%iterator_value)
              Deallocate(comp_buffer)
            End Select
            Call A%ks_points(A%iterator_value)%data%free()
            this_ks = A%iterator_next()
          Enddo
          A%status = STATUS_DISK
        End Select
      !! **********************************************************
      Case(MEMORY_EXECUTION_DISK)
      !! **********************************************************
        If(A%status.ne.STATUS_DISK) Call dash_error(0,"ks_array_to_store","data not on disk")
        Selectcase(memory_storage_mode)
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_DEVICE)
          If(.not.A%storage%initialized()) Call dash_error(0,"ks_array_to_store","disk not initialized")
          this_ks = A%iterator_init()
          Do While(this_ks%exist())
            Call A%ks_points(A%iterator_value)%data%store(1,this_ks%ks_size)
            Select Case(this_ks%ks_type)
            Case(POINT_IS_REAL)
              Allocate(real_buffer(this_ks%ks_size,this_ks%ks_size))
              Call A%storage%data_from_disk(real_buffer,A%iterator_value)
              Call A%ks_points(A%iterator_value)%data%set_raw(real_buffer)
              Deallocate(real_buffer)
            Case(POINT_IS_COMPLEX)
              Allocate(comp_buffer(this_ks%ks_size,this_ks%ks_size))
              Call A%storage%data_from_disk(comp_buffer,A%iterator_value)
              Call A%ks_points(A%iterator_value)%data%set_raw(comp_buffer)
              Deallocate(comp_buffer)
            End Select
            this_ks = A%iterator_next()
          Enddo
          A%status = STATUS_DEVICE
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_HOST)
          If(.not.A%storage%initialized()) Call dash_error(0,"ks_array_to_store","disk not initialized")
          this_ks = A%iterator_init()
          Do While(this_ks%exist())
            Call A%ks_points(A%iterator_value)%data%host_store(1,this_ks%ks_size)
            Select Case(this_ks%ks_type)
            Case(POINT_IS_REAL)
              Allocate(real_buffer(this_ks%ks_size,this_ks%ks_size))
              Call A%storage%data_from_disk(real_buffer,A%iterator_value)
              Call A%ks_points(A%iterator_value)%data%host_set_raw(real_buffer)
              Deallocate(real_buffer)
            Case(POINT_IS_COMPLEX)
              Allocate(comp_buffer(this_ks%ks_size,this_ks%ks_size))
              Call A%storage%data_from_disk(comp_buffer,A%iterator_value)
              Call A%ks_points(A%iterator_value)%data%host_set_raw(comp_buffer)
              Deallocate(comp_buffer)
            End Select
            this_ks = A%iterator_next()
          Enddo
          A%status = STATUS_HOST
        End Select
      !! **********************************************************
      Case(MEMORY_EXECUTION_HOST)
      !! **********************************************************
        If(A%status.ne.STATUS_HOST) Call dash_error(0,"ks_array_to_store","data not on host")
        Selectcase(memory_storage_mode)
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_DEVICE)
          this_ks = A%iterator_init()
          Do While(this_ks%exist())
            Call A%ks_points(A%iterator_value)%data%store(1,this_ks%ks_size)
            Call A%ks_points(A%iterator_value)%data%host_to_dev()
            Call A%ks_points(A%iterator_value)%data%host_free()
            this_ks = A%iterator_next()
          Enddo
          A%status = STATUS_DEVICE
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_DISK)
          If(.not.A%storage%initialized()) Then
            Allocate(all_metatdata(A%n_ks_points))
            this_ks = A%iterator_init()
            Do While (this_ks%exist())
              Call all_metatdata(A%iterator_value)%create( A%iterator_value,A%n_ks_points,this_ks%ks_type,this_ks%ks_size )
              this_ks = A%iterator_next()
            Enddo
            Call A%storage%init('tmp',all_metatdata)
            Deallocate(all_metatdata)
          Endif
          this_ks = A%iterator_init()
          Do While (this_ks%exist())
            Select Case(this_ks%ks_type)
            Case(POINT_IS_REAL)
              Call A%ks_points(A%iterator_value)%data%host_get_raw(real_buffer)
              Call A%storage%data_to_disk(real_buffer,A%iterator_value)
              Deallocate(real_buffer)
            Case(POINT_IS_COMPLEX)
              Call A%ks_points(A%iterator_value)%data%host_get_raw(comp_buffer)
              Call A%storage%data_to_disk(comp_buffer,A%iterator_value)
              Deallocate(comp_buffer)
            End Select
            Call A%ks_points(A%iterator_value)%data%host_free()
            this_ks = A%iterator_next()
          Enddo
          A%status = STATUS_DISK
        End Select
      !! **********************************************************
      End Select
    End Subroutine ks_array_to_store

    Subroutine ks_array_from_store(A)
    !! ============================================================
    !! Procedure used to recover the ks_array form "away" storage
    !! Useful when need to free space on device/host for other
    !! calculations (maybe integrals?)
    !! 
    !! Usage: after the array has been cretated and %global_store-d
    !!        Call A%to_store() will move the array data to storage
    !!        location
    !!        Before using the array, call A%from_store()
    !!
    !! This procedure is general for the whole array, to handle 
    !! single ks points use %load_data/%dump_data
    !!
    !! This is ment to be used after the LA section in the SCF 
    !! (fdik + fermi + pdig)
    !! ============================================================
    Use dash_utils
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Type(info_ks_point)                              :: this_ks
      Real(fp64)         , Dimension(:,:), Allocatable :: real_buffer
      Complex(fp64)      , Dimension(:,:), Allocatable :: comp_buffer
      If(A%size().eq.0) Return
      If(.not.A%initialized()) Call dash_error(0,"ks_array_from_store","array not initialized")
      !! If execution is equal to store nothing to do 
      If(memory_execution_mode.eq.memory_storage_mode) Return

      Select Case(memory_execution_mode)
      !! **********************************************************
      Case(MEMORY_EXECUTION_DEVICE)
      !! **********************************************************
      !! ----------------------------------------------------------
        Selectcase(memory_storage_mode)
        Case(MEMORY_STORAGE_HOST)
          If(A%status.ne.STATUS_HOST) Call dash_error(0,"ks_array_from_store","data not on host")
          this_ks = A%iterator_init()
          Do While(this_ks%exist())
            Call A%ks_points(A%iterator_value)%data%store(1,this_ks%ks_size)
            Call A%ks_points(A%iterator_value)%data%host_to_dev()
            Call A%ks_points(A%iterator_value)%data%host_free()
            this_ks = A%iterator_next()
          Enddo
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_DISK)
          If(A%status.ne.STATUS_DISK) Call dash_error(0,"ks_array_from_store","data not on disk")
          this_ks = A%iterator_init()
          Do While (this_ks%exist())
            Call A%ks_points(A%iterator_value)%data%store(1,this_ks%ks_size)
            Select Case(this_ks%ks_type)
            Case(POINT_IS_REAL)
              Allocate(real_buffer(this_ks%ks_size,this_ks%ks_size))
              Call A%storage%data_from_disk(real_buffer,A%iterator_value)
              Call A%ks_points(A%iterator_value)%data%set_raw(real_buffer)
              Deallocate(real_buffer)
            Case(POINT_IS_COMPLEX)
              Allocate(comp_buffer(this_ks%ks_size,this_ks%ks_size))
              Call A%storage%data_from_disk(comp_buffer,A%iterator_value)
              Call A%ks_points(A%iterator_value)%data%set_raw(comp_buffer)
              Deallocate(comp_buffer)
            End Select
            this_ks = A%iterator_next()
          Enddo
        End Select
        A%status = STATUS_DEVICE
      !! **********************************************************
      Case(MEMORY_EXECUTION_DISK)
      !! **********************************************************
        Selectcase(memory_storage_mode)
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_DEVICE)
          If(A%status.ne.STATUS_DEVICE) Call dash_error(0,"ks_array_from_store","data not on device")
          this_ks = A%iterator_init()
          Do While(this_ks%exist())
            Select Case(this_ks%ks_type)
            Case(POINT_IS_REAL)
              Call A%ks_points(A%iterator_value)%data%get_raw(real_buffer)
              Call A%storage%data_to_disk(real_buffer,A%iterator_value)
              Deallocate(real_buffer)
            Case(POINT_IS_COMPLEX)
              Call A%ks_points(A%iterator_value)%data%get_raw(comp_buffer)
              Call A%storage%data_to_disk(comp_buffer,A%iterator_value)
              Deallocate(comp_buffer)
            End Select
            Call A%ks_points(A%iterator_value)%data%free()
            this_ks = A%iterator_next()
          Enddo
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_HOST)
          this_ks = A%iterator_init()
          Do While(this_ks%exist())
            Select Case(this_ks%ks_type)
            Case(POINT_IS_REAL)
              Call A%ks_points(A%iterator_value)%data%host_get_raw(real_buffer)
              Call A%storage%data_to_disk(real_buffer,A%iterator_value)
              Deallocate(real_buffer)
            Case(POINT_IS_COMPLEX)
              Call A%ks_points(A%iterator_value)%data%host_get_raw(comp_buffer)
              Call A%storage%data_to_disk(comp_buffer,A%iterator_value)
              Deallocate(comp_buffer)
            End Select
            Call A%ks_points(A%iterator_value)%data%host_free()
            this_ks = A%iterator_next()
          Enddo
        End Select
        A%status = STATUS_DISK
      !! **********************************************************
      Case(MEMORY_EXECUTION_HOST)
      !! **********************************************************
        Selectcase(memory_storage_mode)
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_DEVICE)
          If(A%status.ne.STATUS_DEVICE) Call dash_error(0,"ks_array_from_store","data not on device")
          this_ks = A%iterator_init()
          Do While(this_ks%exist())
            Call A%ks_points(A%iterator_value)%data%host_store(1,this_ks%ks_size)
            Call A%ks_points(A%iterator_value)%data%dev_to_host()
            Call A%ks_points(A%iterator_value)%data%free()
            this_ks = A%iterator_next()
          Enddo
      !! ----------------------------------------------------------
        Case(MEMORY_STORAGE_DISK)
          this_ks = A%iterator_init()
          Do While (this_ks%exist())
            Call A%ks_points(A%iterator_value)%data%host_store(1,this_ks%ks_size)
            Select Case(this_ks%ks_type)
            Case(POINT_IS_REAL)
              Allocate(real_buffer(this_ks%ks_size,this_ks%ks_size))
              Call A%storage%data_from_disk(real_buffer,A%iterator_value)
              Call A%ks_points(A%iterator_value)%data%host_set_raw(real_buffer)
              Deallocate(real_buffer)
            Case(POINT_IS_COMPLEX)
              Allocate(comp_buffer(this_ks%ks_size,this_ks%ks_size))
              Call A%storage%data_from_disk(comp_buffer,A%iterator_value)
              Call A%ks_points(A%iterator_value)%data%host_set_raw(comp_buffer)
              Deallocate(comp_buffer)
            End Select
            this_ks = A%iterator_next()
          Enddo
        End Select
        A%status = STATUS_HOST
      !! **********************************************************
      End Select
    End Subroutine ks_array_from_store





    Subroutine ks_array_set_indexing_from_array(A,raw_indeces)
    !! Set or update my_indeces in ks_array
    !! If indeces are already assigned then add new ones after
    Use dash_utils
      Class(ks_array),               Intent(InOut) :: A
      Integer        , Dimension(:), Intent(In   ) :: raw_indeces
      ! Local Variables
      Integer                            :: i
      Integer                            :: idx
      Integer                            :: len
      Integer                            :: old_len = -9
      Integer, Dimension(:), Allocatable :: buffer
      If(raw_indeces(1).eq.POINT_NOT_EXIST) Return
      len = Size(raw_indeces)
      If(len.lt.1) Call dash_error(0,'ks_array_set_indexing','invalid lenght')
      If(Allocated(A%my_indeces)) Then
        old_len = Size(A%my_indeces)
        len = old_len + len
        buffer = A%my_indeces
        Deallocate(A%my_indeces)
      Endif
      Allocate( A%my_indeces( len ) )
      i = 1
      If(Allocated(buffer)) Then
        A%my_indeces(1:old_len) = buffer
        i = i + old_len
      Endif
      Do idx = 1, Size(raw_indeces)
        A%my_indeces(i) = raw_indeces(idx)
        i = i + 1
      Enddo
      A%status = STATUS_INVALID
    End Subroutine ks_array_set_indexing_from_array

    Subroutine ks_array_print_distribution(A,weights,dist_aglo,out)
    !! ============================================================
    !! Print the ks point distribution in out unit
    !! A should be initialized with create/create_vary 
    !! (just %my_indeces and %all_points_info are required)
    !! ============================================================
    Use interface_MPI!, Only: world
    Use backend_module, Only: bk_get_current_device
      Class(ks_array),               Intent(In)           :: A
      Real(fp64)     , Dimension(2), Intent(In), Optional :: weights
      Integer        ,               Intent(In), Optional :: dist_aglo
      Integer        ,               Intent(In), Optional :: out
      ! Local variables
      Integer                                 :: my_rank
      Integer                                 :: i
      Integer                                 :: proc, max_proc
      Integer                                 :: ks, my_ks, max_n_ks
      Real(fp64)                              :: current_weight, max_weight, min_weight
      Integer   , Dimension(:)  , Allocatable :: my_ks_points
      Integer   , Dimension(:)  , Allocatable :: buffer
      Integer   , Dimension(:,:), Allocatable :: ks_data
      Integer   , Dimension(:)  , Allocatable :: buffer_devs
      Integer                                 :: unit = 6
      Integer                                 :: batches
      If(Present(out)) unit = out
      Call world%get(rank=my_rank,size=max_proc)
      If(Allocated(A%my_indeces)) Then
        my_ks = Size(A%my_indeces)
      Else
        my_ks = 0
      Endif
      max_n_ks = my_ks
      Call world%max(max_n_ks)
      Allocate(my_ks_points(max_n_ks))
      my_ks_points = INVALID
      If(Allocated(A%my_indeces)) Then
        my_ks_points(1:my_ks) = A%my_indeces(1:my_ks)
      Endif
      Call world%gather(my_ks_points,buffer)
      Call world%gather(bk_get_current_device(),buffer_devs)
      If(my_rank.eq.0) Then
        Allocate(ks_data(max_n_ks,max_proc))
        ks_data = Reshape(buffer,[max_n_ks,max_proc])
        Write(unit,"(1X,79('*'))")
        Write(unit,"(22X,A)", Advance='no') '*** KS POINTS'
        If(Present(dist_aglo)) Then
          Select Case(dist_aglo)
          Case(KS_DISTRIB_ALG_STATIC)
            Write(unit,"(1X,A)") 'NEW DISTRIBUTION ***'
            If(Present(weights)) Then
              batches = Count( Any(ks_data.ge.0, Dim=1) )
              current_weight = 0.0_fp64
              Do i = 1, Size(A%all_points_info)
                current_weight = current_weight + Merge(weights(1),weights(2),A%all_points_info(i)%ks_type.eq.POINT_IS_REAL)
              Enddo
              Write(unit,"(1X,'KS WEIGHTS:   (R)',F6.2,'   (C)',F6.2,T55,'TARGET WEIGHT:',F12.2)") weights(1:2), current_weight/batches
            Endif
            Write(unit,"(1X,79('*'))")
          Case(KS_DISTRIB_ALG_COMPATIB)
            Write(unit,"(1X,A)") 'LECAGY DISTRIBUTION ***'
            Write(unit,"(1X,79('*'))")
          Case(KS_DISTRIB_ALG_GREED)
            Write(unit,"(1X,A)") 'GREED DISTRIBUTION ***'
            If(Present(weights)) Write(unit,"(1X,'KS WEIGHTS:   (R)',F6.2,'   (C)',F6.2)") weights(1:2)
            Write(unit,"(1X,79('*'))")
          End Select
        Else
          Write(unit,"(1X,A)") 'DISTRIBUTION ***'
          Write(unit,"(1X,79('*'))")
        Endif
        max_weight = 0.0_fp64
        min_weight = Huge(min_weight)
        Do proc = 1, max_proc
          If(ks_data(1,proc).eq.INVALID) Cycle
          current_weight = 0.0_fp64
          i = 1
          Write(unit,"(1X,'P',I5,': ')", Advance='no') proc - 1
          Do ks = 1, max_n_ks
            If(i.gt.8) Then
              Write(unit,*)
              Write(unit,"(1X,A8)",Advance='no') '      : '
              i = 1
            Endif
            If(ks_data(ks,proc).eq.INVALID) Exit
            Associate( this_ks => A%all_points_info(ks_data(ks,proc)))
              Select Case(this_ks%ks_type)
              Case(POINT_IS_REAL)
                Write(unit,"(A1,'(',I4,')')",Advance='no') 'R',this_ks%k_index
                If(Present(weights)) current_weight = current_weight + weights(1)
              Case(POINT_IS_COMPLEX)
                Write(unit,"(A1,'(',I4,')')",Advance='no') 'C',this_ks%k_index
                If(Present(weights)) current_weight = current_weight + weights(2)
              End Select
              Select Case(this_ks%s_index)
              Case(SPIN_ALPHA)
                Write(unit,"(A1,1X)", Advance='no') 'a'
              Case(SPIN_BETA)
                Write(unit,"(A1,1X)", Advance='no') 'b'
              End Select
            End Associate
            i = i + 1
          Enddo
          Write(unit,*)
          If(Present(weights)) Then
            Write(unit,"(1X,A15,F12.2)") '      : Weight ', current_weight
            max_weight = Max(max_weight,current_weight)
            min_weight = Min(min_weight,current_weight)
          Endif
          Write(unit,"(1X,A15,2X,I2)") '      : Device ', buffer_devs(proc)
        Enddo
        If(Present(weights)) Write(unit,"(1X,A,F12.2)") 'Maximum weight difference: ', max_weight - min_weight
        Deallocate(buffer)
        Deallocate(buffer_devs)
        Write(unit,"(1X,79('*'))")
      Endif
    End Subroutine ks_array_print_distribution

    Subroutine ks_array_create(A,source)
    !! ============================================================
    !! Initialize the ks_array from a source
    !! Allocate requred space for the array of points ans sets 
    !! all the info in in ks_point_data
    !!
    !! Do not allocate actual data for matrices
    !! After create need to call global_store to reserve memory
    !! ============================================================
    Use dash_utils
    use interface_MPI
      Class(ks_array), Intent(InOut) :: A
      Type(ks_array) , Intent(In   ) :: source
      ! Local variables
      Logical :: is_complex
      Integer :: nks
      Integer :: ks
      Integer :: i
      Integer :: ierr
      nks = source%n_ks_points      !! points in the mold array
      A%n_ks_points      =  nks
      A%n_ks_points_type =  source%n_ks_points_type
      A%all_points_info  => source%all_points_info
      !! If the array is already created exit
      !! This is necessary as we may call create also for array already created, inside SCF
      If( A%initialized() ) Call dash_error(0,"ks_array_create","array already initialized")
      !! If the array does not contain any ks then return
      If(nks.eq.0) Return
      A%my_indeces = source%my_indeces !! Automatic Fortran reallocation on assignment
      !! Juast a Check
      If(Size(A%my_indeces).ne.nks) Call dash_error(0,'ks_array_create','wrong indeces count',Size(A%my_indeces))
      !! Allocate the required space for ks points
      Allocate(A%ks_points(nks), stat=ierr)
      If(ierr.ne.0) Call dash_error(0,'ks_array_create','allocation A%ks_points failed')
      Loop_on_owned_ks: Do ks = 1, nks
        !! Copy info from all to specific point_data
        A%ks_points(ks)%info => A%all_points_info(A%my_indeces(ks))
        !! Allocate the ks_point (ks_point_data%data) inside each ks_point_data
        is_complex = A%ks_points(ks)%info%ks_type.eq.POINT_IS_COMPLEX
        !! 1 is the number of irreps (NO SYMM: always 1...)
        Call A%ks_points(ks)%data%create(is_complex,1)
        A%ks_points(ks)%status = POINT_STATUS_INIT
        !! At this point the matrix is allocated according to the type of the k point
        !! (real or complex) but the actual data structure is not allocated.
      Enddo Loop_on_owned_ks
      A%status = STATUS_INIT
    End Subroutine ks_array_create

    Subroutine ks_array_create_vary(A,info_array)
    !! ============================================================
    !! Generate a ks_array from a info_ks_point array based on the 
    !! A%my_indeces already assigned
    !!
    !! Just counts how many points are assigned
    !!
    !! >>> KS INDECES OF THE ARRAY SHOULD BE DEFINED BEFOREHAND
    !!
    !! This procedure do not allocate the points vector!
    !! This procedure is used to create the base mold for 
    !! all the other arrays
    !! ============================================================
    Use dash_utils
      Class(ks_array)    ,                       Intent(InOut) :: A
      Type(info_ks_point), Dimension(:), Target, Intent(In   ) :: info_array
      ! Local Variables
      Logical :: is_complex
      Integer :: nks
      Integer :: k, i
      If(.Not.Allocated(A%my_indeces)) Then
        A%all_points_info => info_array
        A%n_ks_points = 0
        A%n_ks_points_type(1) = 0
        A%n_ks_points_type(2) = 0
        Return
      Endif
      nks = Size(A%my_indeces)
      A%all_points_info => info_array
      !! Counting ks types
      Do i = 1, nks
        Select Case(info_array(A%my_indeces(i))%get_type())
        Case(POINT_IS_REAL)
          A%n_ks_points_type(1) = A%n_ks_points_type(1) + 1
          A%n_ks_points = A%n_ks_points + 1
        Case(POINT_IS_COMPLEX)
          A%n_ks_points_type(2) = A%n_ks_points_type(2) + 1
          A%n_ks_points = A%n_ks_points + 1
        Case(POINT_NOT_EXIST)
          !! Should never be here....
          A%n_ks_points = 0
          A%n_ks_points_type(1) = 0
          A%n_ks_points_type(2) = 0
        End Select
      Enddo
      A%status = STATUS_INVALID
    End Subroutine ks_array_create_vary

    Subroutine ks_array_destroy(A)
    !! ============================================================
    !! Deallocate array and destroy
    !! ============================================================
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Type(info_ks_point) :: this_ks
      Call A%global_free()
      A%all_points_info => Null()
      If(Allocated(A%ks_points)) Deallocate(A%ks_points)
      A%n_ks_points = 0
      A%n_ks_points_type = (0,0)
      If(Allocated(A%my_indeces)) Deallocate(A%my_indeces)
      Call A%storage%close()
      A%status = STATUS_INVALID
    End Subroutine ks_array_destroy



!! --- LA Operations ---
    Subroutine ks_array_initialize_ortogonalization_matrices(A)
    !! ============================================================
    !! Construction of first set of ortogonalization matrices
    !! ============================================================
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        !! ********************************************************
        Call A%load_data()
        !! ********************************************************
        Associate(iv => A%iterator_value)
          Call A%ks_points(iv)%data%cholesky()    !! get U from A = U^T * U
          Call A%ks_points(iv)%data%clean_lower() !! set lower triang. to 0
          Call A%ks_points(iv)%data%invert()      !! Compute inverted matrix
        End Associate
        !! ********************************************************
        Call A%dump_data()
        !! ********************************************************
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_initialize_ortogonalization_matrices

    Subroutine ks_array_similarity(A,B,C)
    !! ============================================================
    !! Perform similarity operation on A array, using B array, C is temp
    !!
    !! TODO: WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WKP
    !!
    !! ============================================================
      Class(ks_array), Intent(InOut) :: A
      Type(ks_array) , Intent(InOut) :: B
      Type(ks_array) , Intent(InOut) :: C
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        !! ********************************************************
        Call A%load_data(A%iterator_value)
        Call B%load_data(A%iterator_value)
        Call C%load_data(A%iterator_value,disposable=.True.)
        !! ********************************************************
        Associate(iv => A%iterator_value)
          If(A%ks_points(iv)%data%allocated()) Then
            Call B%ks_points(iv)%data%dagger(OP_DAGGER)
            Call A%ks_points(iv)%data%dagger(OP_NONE)
            Call C%ks_points(iv)%data%multiply(B%ks_points(iv)%data,A%ks_points(iv)%data)
            Call C%ks_points(iv)%data%dagger(OP_NONE)
            Call B%ks_points(iv)%data%dagger(OP_NONE)
            Call A%ks_points(iv)%data%multiply(C%ks_points(iv)%data,B%ks_points(iv)%data)
          Endif
        End Associate
        !! ********************************************************
        Call A%dump_data(A%iterator_value)
        Call A%dump_data(B%iterator_value,disposable=.True.)
        Call A%dump_data(C%iterator_value,disposable=.True.)
        !! ********************************************************
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_similarity

    Subroutine ks_array_multiply(A,B,opb,C,opc)
    !! ============================================================
    !! Multiplication of 2 ks_arrays
    !! ============================================================
      Class(ks_array), Intent(InOut) :: A
      Type(ks_array) , Intent(InOut) :: B
      Integer        , Intent(In   ) :: opb
      Type(ks_array) , Intent(InOut) :: C
      Integer        , Intent(In   ) :: opc
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        !! ********************************************************
        Call A%load_data(A%iterator_value,disposable=.True.)
        Call B%load_data(A%iterator_value)
        Call C%load_data(A%iterator_value)
        !! ********************************************************
        Associate(iv => A%iterator_value)
          If(A%ks_points(iv)%data%allocated()) Then
            Call B%ks_points(iv)%data%dagger(opb)
            Call C%ks_points(iv)%data%dagger(opc)
            Call A%ks_points(iv)%data%multiply(B%ks_points(iv)%data,C%ks_points(iv)%data)
          Endif
        End Associate
        !! ********************************************************
        Call A%dump_data(A%iterator_value)
        Call B%dump_data(A%iterator_value,disposable=.True.)
        Call C%dump_data(A%iterator_value,disposable=.True.)
        !! ********************************************************
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_multiply

    Subroutine ks_array_diag(A,evals)
    !! ============================================================
    !! Diagonalize A array
    !! Returns eigenvalues in ascendig order in evals and
    !! egenvectors in A (same order)
    !! ============================================================
      Class(ks_array)  , Intent(InOut) :: A
      Type(ks_array_1D), Intent(InOut) :: evals
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        !! ********************************************************
        Call A%load_data()
        !! ********************************************************
        Associate( Aks_dat => A%ks_points(A%iterator_value)%data, &
                   e_dativ => evals%data(A%iterator_value)%data   )
          Call Aks_dat%diag(e_dativ)
        End Associate
        !! ********************************************************
        Call A%dump_data()
        !! ********************************************************
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_diag

    Subroutine ks_array_copy(A,B)
    !! ============================================================
    !! Transfer all data stored from B to A, then deallocate B
    !! A should be deallocated ( no global_store called )
    !! ============================================================
    Use dash_utils
      Class(ks_array), Intent(InOut) :: A
      Type(ks_array) , Intent(InOut) :: B
      ! Local Variables
      Type(info_ks_point) :: this_ks
      If(A%size().eq.0) Return
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DISK)
        If( B%status.ne.STATUS_DISK ) Call dash_error(0,"ks_array_copy","data not on disk")
        !! quick way to move disk data associated from B to A
        Call A%storage%transfer_ownership(B%storage)
        A%status = STATUS_DISK
      Case(MEMORY_EXECUTION_DEVICE)
        If( B%status.ne.STATUS_DEVICE ) Call dash_error(0,"ks_array_copy","data not on device")
        this_ks = A%iterator_init()
        Do While(this_ks%exist())
          Associate(iv => A%iterator_value)
              Call A%ks_points(iv)%data%copy(B%ks_points(iv)%data)
          End Associate
          this_ks = A%iterator_next()
        Enddo
        A%status = STATUS_DEVICE
      Case(MEMORY_EXECUTION_HOST)
        If( B%status.ne.STATUS_HOST ) Call dash_error(0,"ks_array_copy","data not on host")
        this_ks = A%iterator_init()
        Do While(this_ks%exist())
          Associate(iv => A%iterator_value)
              Call A%ks_points(iv)%data%host_copy(B%ks_points(iv)%data)
          End Associate
          this_ks = A%iterator_next()
        Enddo
        A%status = STATUS_HOST
      End Select
      b%status = STATUS_INIT
    End Subroutine ks_array_copy




!! --- Iterative Operations ---
    Subroutine ks_array_shift_diagonal(A,shift)
    !! ============================================================
    !! Shift values on the main diagonal, only occupied ones
    !! ============================================================
      Class(ks_array), Intent(InOut) :: A
      Real(fp64)     , Intent(In   ) :: shift
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        If(this_ks%occ_states.gt.0) Then
          !! ******************************************************
          Call A%load_data()
          !! ******************************************************
          Call A%ks_points(A%iterator_value)%data%shift_diagonal(shift,this_ks%occ_states)
          !! ******************************************************
          Call A%dump_data()
          !! ******************************************************
        Endif
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_shift_diagonal


    Subroutine ks_array_set_raw_real(A,raw_data,ind)
    !! ============================================================
    !! Used to store data in array already allocated
    !! Call from ks loop, ind: is the local index of the ks point
    !! raw_adta should be of type coherent with current ks point
    !!
    !! Different behaviours based on memory execution mode:
    !! Place raw_data where we want to store it diring execution
    !!  ... Device ... Host ... Disk ...
    !! ============================================================
      Class(ks_array),                 Intent(InOut) :: A
      Real(fp64)     , Dimension(:,:), Intent(In   ) :: raw_data
      Integer        ,                 Intent(In   ) :: ind
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DISK)
        !! Just store data on disk... so not disturb the device
        Call A%storage%data_to_disk(raw_data,ind)
      Case(MEMORY_EXECUTION_DEVICE)
        !! H2D
        Call A%ks_points(ind)%data%set_raw(raw_data)
      Case(MEMORY_EXECUTION_HOST)
        Call A%ks_points(ind)%data%host_set_raw(raw_data)
      End Select
    End Subroutine ks_array_set_raw_real

    Subroutine ks_array_set_raw_comp(A,raw_data,ind)
    !! ============================================================
    !! Used to store data in array already allocated
    !! Call from ks loop, ind: is the local index of the ks point
    !! raw_adta should be of type coherent with current ks point
    !!
    !! Different behaviours based on memory execution mode:
    !! Place raw_data where we want to store it diring execution
    !!  ... Device ... Host ... Disk ...
    !! ============================================================
      Class(ks_array),                 Intent(InOut) :: A
      Complex(fp64)  , Dimension(:,:), Intent(In   ) :: raw_data
      Integer        ,                 Intent(In   ) :: ind
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DISK)
        Call A%storage%data_to_disk(raw_data,ind)
      Case(MEMORY_EXECUTION_DEVICE)
        !! H2D
        Call A%ks_points(ind)%data%set_raw(raw_data)
      Case(MEMORY_EXECUTION_HOST)
        Call A%ks_points(ind)%data%host_set_raw(raw_data)
      End Select
    End Subroutine ks_array_set_raw_comp

    Subroutine ks_array_get_raw_real(A,raw_data,ind)
    !! ============================================================
    ! Extract raw_data in ks point at position ind
    !! ============================================================
      Class(ks_array),                              Intent(InOut) :: A
      Real(fp64)     , Dimension(:,:), Allocatable, Intent(InOut) :: raw_data
      Integer        ,                              Intent(In   ) :: ind
      ! Local Variables
      Integer :: n, m
      n = A%ks_points(ind)%info%ks_size
      m = A%ks_points(ind)%info%ks_size
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DISK)
        Allocate(raw_data(n,m))
        Call A%storage%data_from_disk(raw_data,ind)
      Case(MEMORY_EXECUTION_DEVICE)
        Call A%ks_points(ind)%data%get_raw(raw_data)
      Case(MEMORY_EXECUTION_HOST)
        Call A%ks_points(ind)%data%host_get_raw(raw_data)
      End Select
    End Subroutine ks_array_get_raw_real

    Subroutine ks_array_get_raw_vector(A,raw_vect,ind)
    !! ============================================================
    !! Return the vector as used by legacy CRYSTAL:
    !! Real ks point:
    !!   * 1D real vector, same memory layout as 2D real array
    !! Complex ks point: 
    !!   * 1D real vector, size = 2*size of 2D complex array
    !!     same memory layout 
    !!       r11,i11,r21,i21,r31,i31,...,r12,i12,r22,i22,...
    !!
    !! The implementation is not fully memory efficient as
    !! all the ...%get_raw procedures will just allocate a 
    !! temporary buffer befor filling the actual raw_vect.
    !! However, using just raw_vect will be much more complicated
    !! so lets stick to this implementation.
    !! ============================================================
    Use dash_utils, Only: dash_error
    Use, Intrinsic :: iso_c_binding, Only: c_loc, c_f_pointer
      Class(ks_array),                            Intent(InOut)           :: A
      Real(fp64)     , Dimension(:), Allocatable, Intent(  Out), Target   :: raw_vect
      Integer        ,                            Intent(In   )           :: ind
      ! Local Variables
      Integer                                      :: n, m
      Integer                                      :: this_ks_type
      Integer                                      :: vec_size
      Complex(fp64), Dimension(:,:), Allocatable :: buffer_comp
      Real(fp64)   , Dimension(:,:), Allocatable :: buffer_real
      Complex(fp64), Dimension(:,:), Pointer     :: vect_comp_pointer

      !! Extra checks in this case as this is a delicate procedure
      !! for interoperability with legacy CRYSTAL code...
      !! Caution is over 9000!!!!
      If(.not.A%initialized()) Call dash_error(0,'ks_array_get_raw_vector','KS Array not allocataed')
      If(ind.gt.Size(A%ks_points)) Call dash_error(0,'ks_array_get_raw_vector','ind too high')
      this_ks_type = A%ks_points(ind)%info%ks_type
      n = A%ks_points(ind)%info%ks_size
      m = A%ks_points(ind)%info%ks_size
      !! Allocate fp64 size if there is at least one complex point
      !! this is the standard behaviour of CRYSTAL
      vec_size = m*n * Merge(2, 1, Any(A%all_points_info(:)%ks_type.eq.POINT_IS_COMPLEX) )
      Allocate(raw_vect(vec_size))
      raw_vect = 0.0_fp64
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DEVICE)
        Select Case(this_ks_type)
        Case(POINT_IS_REAL)
          Call A%ks_points(ind)%data%get_raw(buffer_real)
          raw_vect(1:n*m) = Reshape(buffer_real, [n*m])
          Deallocate(buffer_real)
        Case(POINT_IS_COMPLEX)
          Call A%ks_points(ind)%data%get_raw(buffer_comp)
          Call c_f_pointer(c_loc(raw_vect), vect_comp_pointer, [n,m] )
          vect_comp_pointer = buffer_comp
          Deallocate(buffer_comp)
        End Select
      Case(MEMORY_EXECUTION_DISK)
        Call A%storage%vect_from_disk(raw_vect,ind)
      Case(MEMORY_EXECUTION_HOST)
        Select Case(this_ks_type)
        Case(POINT_IS_REAL)
          Call A%ks_points(ind)%data%host_get_raw(buffer_real)
          raw_vect(1:n*m) = Reshape(buffer_real, [n*m])
          Deallocate(buffer_real)
        Case(POINT_IS_COMPLEX)
          Call A%ks_points(ind)%data%host_get_raw(buffer_comp)
          Call c_f_pointer(c_loc(raw_vect), vect_comp_pointer, [n,m] )
          vect_comp_pointer = buffer_comp
          Deallocate(buffer_comp)
        End Select
      End Select
    End Subroutine ks_array_get_raw_vector

    Subroutine ks_array_get_raw_comp(A,raw_data,ind)
    !! ============================================================
    ! Extract raw_data in ks point at position ind
    !! ============================================================
      Class(ks_array),                              Intent(InOut) :: A
      Complex(fp64)  , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Integer        ,                              Intent(In   ) :: ind
      ! Local Variables
      Integer :: n, m
      n = A%ks_points(ind)%info%ks_size
      m = A%ks_points(ind)%info%ks_size
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DEVICE)
        Call A%ks_points(ind)%data%get_raw(raw_data)
      Case(MEMORY_EXECUTION_DISK)
        Allocate(raw_data(n,m))
        Call A%storage%data_from_disk(raw_data,ind)
      Case(MEMORY_EXECUTION_HOST)
        Call A%ks_points(ind)%data%host_get_raw(raw_data)
      End Select
    End Subroutine ks_array_get_raw_comp

    Function ks_array_iterator_next(A) Result(info)
      Class(ks_array)    , Intent(InOut) :: A
      Type(info_ks_point)                :: info
      A%iterator_value = A%iterator_value + 1
      If(A%iterator_value.le.Ubound(A%ks_points, Dim=1)) Then
        info = A%ks_points(A%iterator_value)%info
      Else
        A%iterator_value = Ubound(A%ks_points, Dim=1) + 1
        info%ks_type = POINT_NOT_EXIST
      Endif
    End Function ks_array_iterator_next

    Subroutine ks_array_iterator_reset(A)
      Class(ks_array), Intent(InOut) :: A
      A%iterator_value = INVALID
    End Subroutine ks_array_iterator_reset

    Function ks_array_iterator_init(A) Result(info)
      Class(ks_array)    , Intent(InOut) :: A
      Type(info_ks_point)                :: info
      A%iterator_value = 1
      If(.not.A%initialized()) Then
        info%ks_type = POINT_NOT_EXIST
        A%iterator_value = INVALID
      Else
        info = A%ks_points(A%iterator_value)%info
      Endif
    End Function ks_array_iterator_init

    !Subroutine ks_array_perform_eigshloc(A,limit,indeces,n_ao)
    !!! This come from inf93, but the options ar not all documented...
    !!! for now not implemented fully....
    !  Class(ks_array),                 Intent(InOut) :: A
    !  Integer        ,                 Intent(In   ) :: limit
    !  Integer        , Dimension(:)  , Intent(In   ) :: indeces
    !  Integer        ,                 Intent(In   ) :: n_ao
    !  Real(fp64)   , Dimension(:,:), Intent(In   ) :: shift
    !  ! Local Variables
    !  Type(info_ks_point) :: this_ks
    !  Integer             :: mu
    !  Integer             :: element
    !  If(l0imit.eq.0) Return
    !  this_ks = A%iterator_init()
    !  Do While(this_ks%exist())
    !    Do mu = 1, limit
    !      element = ((indeces(mu)-1)/n_ao) + 1
    !         Fock( element, element ) = Fock( local_nu, local_mu ) + &
    !              shift(this_ks%s_index,mu)
    !    Enddo
    !    this_ks = A%iterator_next()
    !  Enddo
    !End Subroutine ks_array_perform_eigshloc

    Function ks_array_compute_tester(A) Result (tester)
    !! ============================================================
    !! Compute the tester value, used to check SCF convergence
    !! ============================================================
      Real(fp64)                     :: tester
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Real(fp64)          :: tmp_tester
      Type(info_ks_point) :: this_ks
      tester = 0.0_fp64
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        If(this_ks%occ_states.gt.0) Then
          !! ******************************************************
          Call A%load_data()
          !! ******************************************************
          tmp_tester = A%ks_points(A%iterator_value)%data%compute_tester(this_ks%occ_states)
          tester = Max(tmp_tester,tester)
          !! ******************************************************
          Call A%dump_data()
          !! ******************************************************
        Endif
        this_ks = A%iterator_next()
      Enddo
    End Function ks_array_compute_tester

    Subroutine ks_array_load_data(A,iterator,disposable)
    !! ============================================================
    !! Termporary store ks matrix on device memory, used to allow
    !! different memory execution modes
    !!
    !!           --- USE ONLY INSIDE KS LOOPS ---
    !!
    !! iterator : If not specified, use the current iterator value
    !! disposable : If True DO NOT copy data, just allocate
    !!
    !! * Allocate space on device (%store)
    !! * If data is disposable exit
    !! * Retrive data from disk or host in R/C buffer
    !! * Move buffer to device
    !! ============================================================
    Use dash_utils, Only: dash_error
      Class(ks_array), Intent(InOut)           :: A
      Integer        , Intent(In   ), Optional :: iterator
      Logical        , Intent(In   ), Optional :: disposable
      ! Local Variables
      Real(fp64)   , Dimension(:,:), Allocatable :: real_buffer
      Complex(fp64), Dimension(:,:), Allocatable :: comp_buffer
      Integer                                    :: record_val
      Type(info_ks_point)                        :: this_ks
      If(A%size().eq.0) Return
      If(memory_execution_mode.eq.MEMORY_EXECUTION_DEVICE) Then
        If( A%status.ne.STATUS_DEVICE ) Call dash_error(0,"ks_array_load_data","data not on device")
        Return
      Endif
      record_val = A%iterator_value
      If(Present(iterator)) record_val = iterator
      If(record_val.eq.INVALID) Call dash_error(0,'ks_array_load_data','invalid iterator value')
      this_ks = A%ks_points(record_val)%info
      Call A%ks_points(record_val)%data%store( 1, this_ks%ks_size )
      If(Present(disposable)) Then
        If(disposable) Return
      Endif
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DISK)
        If( A%status.ne.STATUS_DISK ) Call dash_error(0,"ks_array_load_data","data not on disk")
        Select Case(this_ks%ks_type)
        Case(POINT_IS_REAL)
          Allocate(real_buffer(this_ks%ks_size,this_ks%ks_size))
          Call A%storage%data_from_disk(real_buffer,record_val)
          Call A%ks_points(record_val)%data%set_raw(real_buffer)
          Deallocate(real_buffer)
        Case(POINT_IS_COMPLEX)
          Allocate(comp_buffer(this_ks%ks_size,this_ks%ks_size))
          Call A%storage%data_from_disk(comp_buffer,record_val)
          Call A%ks_points(record_val)%data%set_raw(comp_buffer)
          Deallocate(comp_buffer)
        End Select
      Case(MEMORY_EXECUTION_HOST)
        If( A%status.ne.STATUS_HOST ) Call dash_error(0,"ks_array_load_data","data not on host")
        Call A%ks_points(record_val)%data%host_to_dev()
      End Select
    End Subroutine ks_array_load_data

    Subroutine ks_array_dump_data(A,iterator,disposable)
    !! ============================================================
    !! Free termporary storage of ks matrix on device memory, 
    !! used to allow different memory execution modes
    !!
    !!           --- USE ONLY INSIDE KS LOOPS ---
    !!
    !! iterator : If not specified, use the current iterator value
    !! disposable : If .true. DO NOT save data, just deallocate
    !!
    !! * If data is disposable deallocate device memory and exit
    !! * Retrive data from device
    !! * Copy R/C buffer to disk/host
    !! * Deallocate space on device (%free)
    !! ============================================================
    Use dash_utils, Only: dash_error
      Class(ks_array), Intent(InOut)           :: A
      Integer        , Intent(In   ), Optional :: iterator
      Logical        , Intent(In   ), Optional :: disposable
      ! Local Variables
      Real(fp64)   , Dimension(:,:), Allocatable :: real_buffer
      Complex(fp64), Dimension(:,:), Allocatable :: comp_buffer
      Integer                                    :: record_val
      Type(info_ks_point)                        :: this_ks
      If(A%size().eq.0) Return
      If(memory_execution_mode.eq.MEMORY_EXECUTION_DEVICE) Return
      record_val = A%iterator_value
      If(Present(iterator)) record_val = iterator
      If(record_val.eq.INVALID) Call dash_error(0,'ks_array_dump_data','invalid iterator value')
      this_ks = A%ks_points(record_val)%info
      If(Present(disposable)) Then
        If(disposable) Then
          Call A%ks_points(record_val)%data%free()
          Return
        Endif 
      Endif
      Select Case(memory_execution_mode)
      Case(MEMORY_EXECUTION_DISK)
        Select Case(this_ks%ks_type)
        Case(POINT_IS_REAL)
          Call A%ks_points(record_val)%data%get_raw(real_buffer)
          Call A%storage%data_to_disk(real_buffer,record_val)
          Deallocate(real_buffer)
        Case(POINT_IS_COMPLEX)
          Call A%ks_points(record_val)%data%get_raw(comp_buffer)
          Call A%storage%data_to_disk(comp_buffer,record_val)
          Deallocate(comp_buffer)
        End Select
      Case(MEMORY_EXECUTION_HOST)
        Call A%ks_points(record_val)%data%dev_to_host()
      End Select
      Call A%ks_points(record_val)%data%free()
    End Subroutine ks_array_dump_data
















    Subroutine ks_array_compute_density(A,pg_size,density_matrix,force_calc,A1D,shift)
    !! ============================================================
    !! Compute density matrix from eigenvectors array
    !! Eigenvectors should be alredy computed and stored (A)
    !!
    !! For SCF:             force_calc = False
    !!     - compute alpha/beta spins separately, then merge
    !!       as alpha+beta alpha-beta
    !!
    !! For forces:          force_calc = True
    !!     - compute only alpha+beta
    !!     - also weight couple contribution for eigenvalues
    !!
    !! ONLY WORKS ON SYMMREMO, close/open shell, real/complex points
    !!
    !! This procedure automatically allocates the density matrix 
    !! (pg_irr) of the correct size (as CRYSTAL)
    !! ============================================================
    Use dash_utils   , Only: dash_utils_timer,dash_error
    Use interface_MPI, Only: world
      Class(ks_array)  ,                            Intent(InOut)           :: A
      Integer          ,                            Intent(In   )           :: pg_size
      Real(fp64)       , Dimension(:), Allocatable, Intent(InOut), Target   :: density_matrix
      Logical          ,                            Intent(In   )           :: force_calc
      Type(ks_array_1D),                            Intent(InOut), Optional :: A1D
      Real(fp64)       ,                            Intent(In   ), Optional :: shift
      ! Local Variables
      Type(info_ks_point)                           :: this_ks
      Integer                                       :: occupied_states
      Integer                                       :: spin_states
      Logical                                       :: do_alpha, do_beta
      Real(fp64), Dimension(:), Allocatable         :: weights
      Real(fp64), Dimension(:), Allocatable         :: evals
      Real(fp64)                                    :: unit_occupation  !! unitary occupation of bands
      Real(fp64), Dimension(:), Allocatable, Target :: p_alpha
      Real(fp64), Dimension(:), Allocatable, Target :: p_beta
      Real(fp64), Dimension(:), Pointer             :: pg
      spin_states  = Maxval(A%all_points_info(:)%s_index)
      If(force_calc) Then
          Allocate(density_matrix(pg_size))
          density_matrix = 0.0_fp64
          If(.not.Present(A1D)) Call dash_error(0,"ks_array_compute_density","1D array required for force calc.")
          If(.not.Present(shift)) Call dash_error(0,"ks_array_compute_density","shift required for force calc.")
      Else
          do_alpha = Count(A%ks_points(:)%info%s_index.eq.SPIN_ALPHA) .gt. 0
          do_beta  = Count(A%ks_points(:)%info%s_index.eq.SPIN_BETA ) .gt. 0
          If(do_alpha) Then
            Allocate(p_alpha(pg_size))
            p_alpha = 0.0_fp64
          Endif
          If(do_beta) Then
            Allocate(p_beta (pg_size))
            p_beta = 0.0_fp64
          Endif
      Endif
      unit_occupation = Merge( 2._fp64, 1._fp64, spin_states.eq.CLOSE_SHELL ) !! warning for SOC
      Call dash_utils_timer('PDIG_pre',3)
      this_ks = A%iterator_init()
      Do While (this_ks%exist())
        If(this_ks%occ_states.eq.0) Cycle !! Should be impossible to have a k not occupied....
        occupied_states = Count( Abs( this_ks%occ_bands ) .ge. 1.0E-11_fp64 )
        Allocate(weights(occupied_states))
        If(force_calc) Then
          Call A1D%data(A%iterator_value)%data%get_raw(evals)
          weights = this_ks%occ_bands( 1:occupied_states ) * unit_occupation * ( evals(1:occupied_states) + shift )
        Else
          weights = this_ks%occ_bands( 1:occupied_states ) * unit_occupation
        Endif
        !! ********************************************************
        Call A%load_data()
        !! ********************************************************
        If(force_calc) Then
          pg => density_matrix
        Else
          Select Case(this_ks%s_index)
          Case(SPIN_ALPHA)
            pg => p_alpha
          Case(SPIN_BETA)
            pg => p_beta
          End Select
        Endif

        Call A%ks_points(A%iterator_value)%data%compute_density( &
             & occupied_states,pg,weights,this_ks%coord )
        !! ********************************************************
        Call A%dump_data()
        !! ********************************************************
        Deallocate(weights)
        this_ks = A%iterator_next()
      Enddo
      Call dash_utils_timer('PDIG_ks',3)
      If(.not.force_calc) Then
        Select Case(spin_states)
        Case(CLOSE_SHELL)
          Allocate(density_matrix(pg_size))
          density_matrix = 0.0_fp64
          If(do_alpha) density_matrix(1:pg_size) = p_alpha(1:pg_size)
        Case(OPEN_SHELL)
          Allocate(density_matrix(pg_size*2))
          density_matrix = 0.0_fp64
          If(do_alpha) Then
            density_matrix(1:pg_size) = p_alpha(1:pg_size)
            density_matrix(pg_size+1:pg_size*2) = p_alpha(1:pg_size)
          Endif
          If(do_beta ) Then
            density_matrix(1:pg_size) = density_matrix(1:pg_size) + p_beta(1:pg_size)
            density_matrix(pg_size+1:pg_size*2) = density_matrix(pg_size+1:pg_size*2) - p_beta(1:pg_size)
          Endif
        End Select
        Call dash_utils_timer('PDIG_ab',3)
      Endif
      Call world%sum(density_matrix)
      Call dash_utils_timer('PDIG_mpi',3)
    End Subroutine ks_array_compute_density
!!#############################################################################





!!#############################################################################
!!        ks_array_1D type-bound procedures
!!#############################################################################
    ! Public implementations
    Subroutine ks_array_1D_create(A1D,A,ldim)
    !! Allocate the space for a 1D array
    !! This array is always stored in memory, as it has order 1
      Class(ks_array_1D), Intent(InOut) :: A1D
      Type(ks_array)    , Intent(InOut) :: A
      Integer           , Intent(In   ) :: ldim
      ! Local Variables
      Integer :: nks
      Type(info_ks_point) :: this_ks
      !! If this 1D array is already allocated return
      If(A1D%status.eq.POINT_STATUS_ALLOC) Return
      nks = A%n_ks_points
      !! Check if there are points in this array
      If(nks.eq.0) Then
        !Allocate( A1D%data(1:1) )
        !Allocate( A1D%data(1)%info )
        !A1D%data(1)%info%ks_type = POINT_NOT_EXIST
        A1D%all_points_info => all_ks_points_info_target
        A1D%status = POINT_STATUS_ALLOC
        Return
      Endif
      A1D%all_points_info => all_ks_points_info_target
      Allocate(A1D%data(1:nks))
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        Associate( iv => A%iterator_value )
          Call A1D%data(iv)%data%create(ldim)
          A1D%data(iv)%info => A%ks_points(iv)%info
        End Associate
        this_ks = A%iterator_next()
      Enddo
      A1D%status = POINT_STATUS_ALLOC
    End Subroutine ks_array_1D_create

    Subroutine ks_array_1D_destroy(A1D)
    !! Deallocate the space for 1D array
      Class(ks_array_1D), Intent(InOut) :: A1D
      ! Local Variables
      Type(info_ks_point) :: this_ks
      If(A1D%status.eq.POINT_STATUS_ALLOC) Then
        this_ks = A1D%iterator_init()
        Do While(this_ks%exist())
          Associate( iv => A1D%iterator_value )
            Call A1D%data(iv)%data%destroy()
          End Associate
          this_ks = A1D%iterator_next()
        Enddo
        If(Allocated(A1D%data)) Deallocate(A1D%data)
        A1D%all_points_info => Null()
        A1D%status = POINT_STATUS_INIT
      Endif
    End Subroutine ks_array_1D_destroy

    Subroutine ks_array_1D_get_raw(A1D,raw_data,ind)
      Class(ks_array_1D),                            Intent(InOut) :: A1D
      Real(fp64)        , Dimension(:), Allocatable, Intent(  Out) :: raw_data
      Integer           ,                            Intent(In   ) :: ind
      Call A1D%data(ind)%data%get_raw(raw_data)
    End Subroutine ks_array_1D_get_raw

    Subroutine ks_array_1D_set_raw(A1D,raw_data,ind)
      Class(ks_array_1D),               Intent(InOut) :: A1D
      Real(fp64)        , Dimension(:), Intent(In   ) :: raw_data
      Integer           ,               Intent(In   ) :: ind
      Call A1D%data(ind)%data%set_raw(raw_data)
    End Subroutine ks_array_1D_set_raw

    Subroutine ks_array_1D_shift(A1D,shift)
      Class(ks_array_1D), Intent(InOut) :: A1D
      Real(fp64)        , Intent(In   ) :: shift
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A1D%iterator_init()
      Do While (this_ks%exist())
        Associate( occ => this_ks%occ_states, iv => A1D%iterator_value )
          If(occ.gt.0) Call A1D%data(iv)%data%shift(occ,shift)
        End Associate
        this_ks = A1D%iterator_next()
      Enddo
      Call A1D%iterator_reset()
    End Subroutine ks_array_1D_shift

    Subroutine ks_array_1D_to_ene(A1D,ene)
    Use interface_MPI, Only: world
    Use dash_utils
      Class(ks_array_1D),               Intent(InOut) :: A1D
      Real(fp64)        , Dimension(:), Intent(InOut) :: ene
      ! Local Variables
      Integer                               :: n_spins
      Integer                               :: idx, start, finish
      Type(info_ks_point)                   :: this_ks
      Real(fp64), Dimension(:), Allocatable :: buffer
      ene = 0.0_fp64
      n_spins = Maxval(A1D%all_points_info(:)%s_index) 
      this_ks = A1D%iterator_init()
      Do While (this_ks%exist())
        Call A1D%data(A1D%iterator_value)%data%get_raw(buffer)
        idx = (this_ks%k_index - 1) * n_spins + this_ks%s_index
        start  = (idx-1) * this_ks%ks_size + 1
        finish =  idx    * this_ks%ks_size
        ene(start:finish) = buffer
        Deallocate(buffer)        
        this_ks = A1D%iterator_next()
      Enddo
      Call A1D%iterator_reset()
      Call world%sum(ene)
    End Subroutine ks_array_1D_to_ene

    Subroutine ks_array_1D_print(A1D,ene,k_prt,out)
      Class(ks_array_1D),               Intent(InOut) :: A1D
      Real(fp64)        , Dimension(:), Intent(InOut) :: ene
      Integer           ,               Intent(In   ) :: k_prt
      Integer           ,               Intent(In   ) :: out
      ! Local Variables
      Integer             :: i, s, k
      Integer             :: n_spins, n_k_points
      Integer             :: idx, start, finish
      Type(info_ks_point) :: this_ks
      n_spins = Maxval(A1D%all_points_info(:)%s_index)
      n_k_points = Maxval(A1D%all_points_info(:)%k_index)
      i = 1
      Do s = 1, n_spins
        Write(out,"(/A24)") spin_str(s)
        Do k = 1, Min(k_prt,n_k_points)
          this_ks = A1D%all_points_info(i)
          idx = (this_ks%k_index - 1) * n_spins + this_ks%s_index
          start  = (idx-1) * this_ks%ks_size + 1
          finish =  idx    * this_ks%ks_size
          Write(out,"(/I4,'(',3I2,')'/(T13,1P,5(E22.14)))") &
          & this_ks%k_index,this_ks%coord(1:3),ene(start:finish)
        Enddo
      Enddo
    End Subroutine ks_array_1D_print




  Subroutine ks_array_1D_compute_fermi(A1D,nelec,energies,is_shifted,shift,do_smear,smear,band_gap,out)
  !! WARNING: WIP WIP WIP WIP WIP WIP WIP WIP 
    Use interface_MPI
    Use backend_module
    Use dash_utils
      Class(ks_array_1D),               Intent(InOut)         :: A1D
      Integer           ,               Intent(In   )         :: nelec
      Real(fp64)        , Dimension(:), Intent(InOut), Target :: energies
      Logical           ,               Intent(In   )         :: is_shifted
      Real(fp64)        ,               Intent(In   )         :: shift
      Logical           ,               Intent(In   )         :: do_smear
      Real(fp64)        ,               Intent(In   )         :: smear
      Real(fp64)        , Dimension(2), Intent(InOut)         :: band_gap
      Integer           ,               Intent(In   )         :: out
      ! Local Variables
      Integer                  :: nk, ns, nb, k, s, band
      Integer                  :: e, e_alpha, e_beta
      Integer                  :: k_top, k_bottom

      Real(fp64)               :: fermi_guess
      Real(fp64), Dimension(2) :: ranges
      Real(fp64)               :: toll
      Real(fp64)               :: homo,lumo
      Type(info_ks_point)      :: this_ks
      Integer   , Dimension(2) :: full_bands
      Integer   , Dimension(2) :: cutting_bands
      Real(fp64), Dimension(:,:,:), Pointer :: e3d => Null()

      Write(out,"(20('='),1X,A15,1X,20('='))") "START FERMI GPU"
      If(A1D%status.ne.POINT_STATUS_ALLOC) Call dash_error(0,"ks_array_1D_compute_fermi","not allocated")

      !! --- Retrieve base parameters ---
      nk = Count(A1D%all_points_info(:)%s_index.eq.SPIN_ALPHA)
      ns = Maxval(A1D%all_points_info(:)%s_index)
      nb = Maxval(A1D%all_points_info(:)%ks_size)
      
      !! --- Initialize and remove parallelism ---
      energies = 0.0_fp64
      Call A1D%to_ene(energies)
      Call world%sum(energies)
      e3d(1:nb, 1:ns, 1:nk) => energies
      toll = Merge(1.0E-9_fp64, 1.0E-5_fp64, nk.eq.1)
      If( (nk.ne.1) .and. do_smear ) toll = Max(24.0_fp64*smear, toll)

      !! --- Compute homo and lumo values ---
      homo = -Huge(homo)
      lumo =  Huge(lumo)
      If(ns.eq.OPEN_SHELL) Then
        Do k = 1, nk
          e_alpha = 0
          e_beta  = 0
          Do e = 1, nelec !! Loop on electrons and assign them to lowest free band
            If( e_alpha .ge. nb ) Then        !! alpha channel is full
              e_beta = e_beta + 1
            Else If( e_beta .ge. nb ) Then    !! Beta channel is full
              e_alpha = e_alpha +1
            !! No channel is full, so take the lowest
            Else If( e3d(e_alpha+1,1,k) .le. e3d(e_beta+1,2,k) ) Then 
              e_alpha = e_alpha + 1
            Else
              e_beta = e_beta + 1
            Endif
          Enddo !! e loop
          !! Now e_alpha and e_beta have the homos
          homo = Max( homo, e3d(e_alpha,1,k), e3d(e_beta,2,k) )

          e_alpha = Min( e_alpha+1, nb )
          e_beta  = Min( e_beta+1 , nb )
          !! Now e_alpha and e_beta have the lumos
          lumo = Min( lumo, e3d(e_alpha,1,k), e3d(e_beta,2,k) )
        Enddo !! k loop
      Else If(ns.eq.CLOSE_SHELL) Then
        e_alpha = nelec/2 + Mod(nelec,2)  !! total alpha electrons
        homo = Maxval( e3d(e_alpha,1,:) )
        band = (nelec/2)+1                   !! first unoccupied (or not fully occupied) band
        lumo = Minval( e3d(band,1,:) )
      Endif

      !! ---  ---
      fermi_guess = homo
      homo = homo + toll
      lumo = lumo - toll

      Write(out,*) "GET DATA ----------------------------------------------"
      Write(out,*) "MAX   : ",lumo
      Write(out,*) "GUESS : ",fermi_guess
      Write(out,*) "MIN   : ",homo
      Write(out,*) "*******************************************************"


      !! --- Band Analysis ---
      full_bands = 0
      cutting_bands = 0
      Do s = 1, ns            !! spin loop
        Do band = 1, nb       !! loop on bnds
          If( Maxval(e3d(band,s,:)) .le. lumo ) Then
            !! --- below LUMO-toll ---
            full_bands(s) = full_bands(s) + 1
          Else
            !! --- above LUMO-toll, but below HOMO+toll (cutting fermi) ---
            If( Minval(e3d(band,s,:)) .lt. homo ) then
              cutting_bands(s) = cutting_bands(s) + 1
            endif
          Endif
        Enddo
      Enddo
      !NCUTT = Sum(cutting_bands)
      !NCUT  = full(s) + cutting(s)
      !NFUL  = full(s)

      ! --- Set Guess Occupation in global memory ---
      Do k = 1, Size( A1D%all_points_info )
          If(A1D%all_points_info(k)%weight .eq. INVALID ) Call dash_error(0,"compute_fermi","weight not avail")
          A1D%all_points_info(k)%occ_bands = 0._fp64
          Associate( k_ => A1D%all_points_info(k)%k_index, &
                   & s_ => A1D%all_points_info(k)%s_index  )
            A1D%all_points_info(k)%occ_bands( 1:full_bands(s_) ) = A1D%all_points_info(k)%weight
          End Associate
      Enddo

      !! --- Compute bandgap ---
      band_gap = 0.0_fp64
      If( Sum(cutting_bands).eq.0 ) Then

        !If(do_smear) ENTROP(:) = 0._float !jb fix nonzero entropy for insulators

      !  if(skipgap)then !jb: do not print gap if spin locking is active, as it's meaningless
      !    write(iout,"(' SPIN LOCKING: NO ENERGY GAP COMPUTED')")
      !    bandgap(:)=1000._float
      !  else

        Do s = 1, ns
          If( ns.eq.OPEN_SHELL ) Write(out,"(/,A24)") spin_str(s)

          If(is_shifted) Then
            Write(out,"(1X,A33,F7.2,' au')") 'INSULATING STATE - LEVEL SHIFTER ', shift
          Else
            Write(out,"(1X,A)") 'INSULATING STATE'
          Endif

          band = full_bands(s)
          k_top = Maxloc( e3d(band,s,:), Dim=1 )

          Write(out,"(1X,A,T33,I6,'; K ',I4,'; EIG',1PE15.7,' AU')") &
               & 'TOP OF VALENCE BANDS -    BAND', band,             &
               & k_top, Maxval( e3d(band,s,:) )

          If(k_top.ne.1) Write(out,"(1X,A,T33,I6,'; K ',I4,'; EIG',1PE15.7,' AU')") &
                                 & 'TOP OF VALENCE BANDS -    BAND', band, 1, e3d(band,s,1)

          band = band + 1
          If(band.gt.nb) Cycle
          k_bottom = Minloc( e3d(band,s,:), Dim=1 )
          Write(out,"(1X,A,T33,I6,'; K ',I4,'; EIG',1PE15.7,' AU')") &
               & 'BOTTOM OF VIRTUAL BANDS - BNAD', band,             &
               & k_bottom, Minval( e3d(band,s,:) )

          If(is_shifted) Then
            band_gap(s) = ( e3d(band,s,k_bottom) - e3d(band-1,s,k_top) - shift ) * ha_to_ev
          Else
            band_gap(s) = ( e3d(band,s,k_bottom) - e3d(band-1,s,k_top)  ) * ha_to_ev
          Endif

          If( k_top.eq.k_bottom ) Then
            Write(out,"(' DIRECT ENERGY BAND GAP:',1X,1F8.4,' eV')") band_gap(s)
          Else
            Write(out,"(' INDIRECT ENERGY BAND GAP:',1X,1F8.4,' eV')") band_gap(s)
          Endif
          If(k_bottom.ne.1) Write(out,"(1X,A,T33,I6,'; K ',I4,'; EIG',1PE15.7,' AU')") &
                                   & 'BOTTOM OF VIRTUAL BANDS - BNAD', band, 1, e3d(band,s,1)
        Enddo
      Else 
        print*,'===> Conducting system'
        !! CONDUCTING

        If(nk.eq.1) Call dash_error(1,'FERMI','ONE K POINT NOT SUFFICIENT FOR FERMI ENERGY CALCULATION')

      !CALL AB          this do not touch ene
      !  CALL ALLOC_FERMI(NCUTT*NKIF,NKIF,NKF,LSWP)
      !  NBAND=NSPSTA*NDF
      !  ISHF=0
      !  DO ISIGMA=1,NSPSTA                 !! spin loop
      !    DO N=NFUL(ISIGMA)+1,NCUT(ISIGMA) !! band loop
      !      NKN=N
      !      DO NK=1,NKF                    !! k loop
      !        EPS(NK)=ENE(NKN)
      !        NKN=NKN+NBAND
      !      ENDDO
      !      CALL OMEGA
      !    ENDDO
      !  ENDDO
      !  IF(PAR(50).GT.EMI.AND.PAR(50).LT.EM)THEN
      !    PAR(7)=PAR(50)
      !  ELSE
      !    PAR(7)=(EM+EMI)*0.5_FLOAT
      !  ENDIF
      !  IF(LPRINT(78).NE.0) WRITE(IOUT,97)EMI,EM
      !  emiro=emi-.5_FLOAT
      !  emro=em+.5_FLOAT
      !  CALL ZERO(EMIro,EMro,PAR(7),RES,NV)
      !  WRITE(IOUT,101)PAR(7),RES,NV
      !  PAR(50)=PAR(7)
      !  CALL CALPES
      Endif

!C... DETERMINE OCCUPATION PATTERNS FOR EACH K-POINT

!      efermi=par(50)
!      N0=0
!      NK=0
!      IF(.NOT.LSOC)THEN
!       DO K=1,NKF
!       DO ISIGMA=1,NSPSTA
!
!       DO N=IBVAL,NDF
!       IF(ALFA(N0+N).LT.(1E-15_FLOAT))GOTO 4002
!       ENDDO
!       N=NDF+1
!4002   NK=NK+1
!       JALPHA(NK)=N-1
!       N0=N0+NDF
!       enddo
!       enddo
!      ELSE
!       NDF=NDF/2
!       NSPSTA=2
!       DO K=1,NKF
!       DO ISIGMA=1,NSPSTA
!       DO N=IBVAL,2*NDF-1,2
!       IF(ALFA(N0+N+ISIGMA-1).LT.(1E-15_FLOAT))GOTO 4003
!       ENDDO
!       N=2*NDF+1
!4003   NK=NK+1
!       JALPHA(NK)=(N-1)/2
!       enddo
!       N0=N0+2*NDF
!       enddo
!      ENDIF
!      IF (LPRINT(78).GT.0) THEN
!      NK=0
!      DO 4004 K=1,NKF
!      WRITE(IOUT,1984)K,(JALPHA(NK+N),N=1,NSPSTA)
!4004  NK=NK+NSPSTA
!      ENDIF
!97    FORMAT(' LIMITS FOR FERMI ENERGY SEARCHING: ',1P,2E11.3)
!101   FORMAT(' POSSIBLY CONDUCTING STATE - EFERMI(AU)',1P,E15.7,' (RES. CHARGE',E10.2,';IT.',I3,')')
!1984  FORMAT(' (by occupation) K,NALPHA,(NBETA)',3I4)

      Write(out,"(20('='),1X,A15,1X,20('='))") " END FERMI GPU "
    End Subroutine ks_array_1D_compute_fermi



    Subroutine ks_array_1D_sort(A1D,A)
      Class(ks_array_1D), Intent(InOut)           :: A1D
      Type(ks_array)    , Intent(InOut), Optional :: A
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A1D%iterator_init()
      Do While(this_ks%exist())
          If(Present(A)) Then
            !! ****************************************************
            Call A%load_data(A1D%iterator_value)
            !! ****************************************************
            Call A1D%data(A1D%iterator_value)%data%sort(A%ks_points(A1D%iterator_value)%data)
            !! ****************************************************
            Call A%dump_data(A1D%iterator_value)
            !! ****************************************************
          Else
            Call A1D%data(A1D%iterator_value)%data%sort()
          Endif
        this_ks = A1D%iterator_next()
      Enddo

    End Subroutine ks_array_1D_sort















    Function ks_array_1D_iterator_next(A1D) Result(info)
      Class(ks_array_1D) , Intent(InOut) :: A1D
      Type(info_ks_point)                :: info
      A1D%iterator_value = A1D%iterator_value + 1
      If(A1D%iterator_value.le.Ubound(A1D%data, Dim=1)) Then
        info = A1D%data(A1D%iterator_value)%info
      Else
        A1D%iterator_value = Ubound(A1D%data, Dim=1) + 1
        info%ks_type = POINT_NOT_EXIST
      Endif
    End Function ks_array_1D_iterator_next

    Subroutine ks_array_1D_iterator_reset(A1D)
      Class(ks_array_1D), Intent(InOut) :: A1D
      A1D%iterator_value = INVALID
    End Subroutine ks_array_1D_iterator_reset

    Function ks_array_1D_iterator_init(A1D) Result(info)
    Use dash_utils
      Class(ks_array_1D) , Intent(InOut) :: A1D
      Type(info_ks_point)                :: info
      If(A1D%status.ne.POINT_STATUS_ALLOC) Call dash_error(0,"ks_array_1D_iterator_init","ARRAY NOT ALLOCATED")
      If(.not.Allocated(A1D%data)) Then
        info%ks_type = POINT_NOT_EXIST
        Return
      Endif
      A1D%iterator_value = 1
      info = A1D%data(1)%info
    End Function ks_array_1D_iterator_init
!!#############################################################################















End Module module_arrays
