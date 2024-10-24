Module module_arrays

!===========================================================
!    @author: Giacomo Ambrogio       date: Feb 2024
!
!    DASH - Device Accelerated Solutions for HPC
!
!    K point array main class
!
!    This module contains the othermost data structures
!    for dealing with k points
!===========================================================

Use base_numbers
Use module_points

  Implicit None
  Private

  Integer, Parameter, Public :: N_IRREPS_NOT_DEFINED = -1
  Integer, Public :: memory_execution_mode = 2              !! Do not change

  Type, Public :: info_ks_point
    Integer,               Private :: k_index
    Integer,               Private :: s_index
    Integer,               Private :: ks_type
    Integer,               Private :: ks_size
    Integer,               Private :: irrep   = INVALID    !! Not yet implemented
    Integer, Dimension(3), Private :: coord
    Integer,               Private :: nocc
    Contains
    Procedure, Public :: create    => info_ks_point_create




    Procedure, Public :: exist     => info_ks_point_exist
    Procedure, Public :: get_type  => info_ks_point_get_type
    Procedure, Public :: get_coord => info_ks_point_get_coord
    Procedure, Public :: get_spin  => info_ks_point_get_spin
    Procedure, Public :: get_k     => info_ks_point_get_k
    Procedure, Public :: get_nocc  => info_ks_point_get_nocc
  End Type

  Type, Public :: ks_point_data
  !! just a wrapper to contain the info and allocation status of the point
  !! Maybe this is redundant...
    Type(info_ks_point), Private :: info
    Integer            , Private :: status = POINT_STATUS_DEFAULT
    Type(ks_point)     , Private :: data
  End Type

  Type, Public :: ks_array_1D
  !! Type to contains 1D arrays, follows the same structure
  !! as the main array class, but is slim
    Type(ks_point_1D)  , Dimension(:), Allocatable, Private :: data
    Integer            ,                            Private :: status = POINT_STATUS_DEFAULT
    Type(info_ks_point), Dimension(:), Allocatable, Private :: info
    Contains
    Procedure, Public :: create  => ks_array_1D_create
    Procedure, Public :: destroy => ks_array_1D_destroy
    Procedure, Public :: get_raw => ks_array_1D_get_raw
    Procedure, Public :: set_raw => ks_array_1D_set_raw
  End Type

  Type, Public :: ks_array
    Type(info_ks_point), Dimension(:), Allocatable, Private :: all_points_info
    Type(ks_point_data), Dimension(:), Allocatable, Private :: ks_points
    Integer                                       , Private :: n_ks_points      = 0
    Integer            , Dimension(2)             , Private :: n_ks_points_type = (0,0)
    Integer            , Dimension(:), Allocatable, Private :: my_indeces
    Integer                                       , Private :: iterator_value   = INVALID
    Contains
      !! --- LA Operations ---
    Procedure, Public :: initialize_ortog_mat => ks_array_initialize_ortogonalization_matrices
    Procedure, Public :: multiply             => ks_array_multiply
    Procedure, Public :: diag                 => ks_array_diag
    Procedure, Public :: copy                 => ks_array_copy
      !! --- Iterative Operations ---
    Procedure, Public :: shift_diagonal => ks_array_shift_diagonal
      !! --- General Operations ---
    Procedure, Public :: set_nocc => ks_array_set_nocc
      !! --- Utilities ---
    Generic  , Public :: create       => ks_array_create
    Generic  , Public :: create       => ks_array_create_vary
    Procedure, Public :: global_store => global_store_ks_array
    Procedure, Public :: global_free  => global_free_ks_array
    Procedure, Public :: allocated    => ks_array_allocated
    Procedure, Public :: initialized  => ks_array_initialized
    Procedure, Public :: size         => ks_array_get_number_of_ks_points
    Procedure, Public :: set_indexing => ks_array_set_indexing_from_array
    Procedure, Public :: print_dist   => ks_array_print_distribution
    Generic  , Public :: set_raw      => set_raw_real, set_raw_comp
    Generic  , Public :: get_raw      => get_raw_real, get_raw_comp, get_raw_vect
    Procedure, Public :: compute_tester          => ks_array_compute_tester
    Procedure, Public :: iterator_init           => ks_array_iterator_init
    Procedure, Public :: iterator_next           => ks_array_iterator_next
    Procedure, Public :: iterator_prev           => ks_array_iterator_prev
    Procedure, Public :: iterator_reset          => ks_array_iterator_reset
    Procedure, Public :: local_store             => local_store_ks_array
    Procedure, Public :: local_free              => local_free_ks_array
    ! Private Implementations
    Procedure, Private :: ks_array_create
    Procedure, Private :: ks_array_create_vary
    Procedure, Private :: move                    => ks_array_move
    Procedure, Private :: set_raw_real            => ks_array_set_raw_real
    Procedure, Private :: set_raw_comp            => ks_array_set_raw_comp
    Procedure, Private :: get_raw_real            => ks_array_get_raw_real
    Procedure, Private :: get_raw_vect            => ks_array_get_raw_vector
    Procedure, Private :: get_raw_comp            => ks_array_get_raw_comp
    Procedure, Private :: initialize_ks_array_from_source
    Procedure, Private :: initialize_ks_array
    Procedure, Private :: count_point                     => increment_ks_points_counter
  End Type

  Contains
!!#############################################################################
!!        info_ks_point type-bound procedures
!!#############################################################################
    Subroutine info_ks_point_create(I,k_index,s_index,ks_type,ks_size,n_of_irreps,n_occ_ao,k_coords)
    !! Reset all info of this ks, and sets new if present
    !! Used to generete the first set of info_ks_point
    Use dash_utils, Only: dash_error
      Class(info_ks_point) , Intent(InOut)           :: I
      Integer              , Intent(In   ), Optional :: k_index
      Integer              , Intent(In   ), Optional :: s_index
      Integer              , Intent(In   ), Optional :: ks_type
      Integer              , Intent(In   ), Optional :: ks_size
      Integer              , Intent(In   ), Optional :: n_of_irreps
      Integer              , Intent(In   ), Optional :: n_occ_ao
      Integer, Dimension(3), Intent(In   ), Optional :: k_coords
      If(Present(n_of_irreps).and.(n_of_irreps.gt.1)) Then
        Call dash_error(0,"info_ks_point_create","number of irreps > 1")
      Endif
      I%k_index = POINT_NOT_EXIST
      I%s_index = SPIN_ALPHA
      I%ks_type = POINT_NOT_EXIST
      I%ks_size = 1
      I%irrep   = 1
      I%coord   = [0,0,0]
      I%nocc    = 0
      If(Present(k_index    )) I%k_index = k_index
      If(Present(s_index    )) I%s_index = s_index
      If(Present(ks_type    )) I%ks_type = ks_type
      If(Present(ks_size    )) I%ks_size = ks_size
      If(Present(n_of_irreps)) I%irrep   = n_of_irreps
      If(Present(k_coords   )) I%coord   = k_coords
      If(Present(n_occ_ao   )) I%nocc    = n_occ_ao
    End Subroutine info_ks_point_create

!!#############################################################################
!!        ks_array type-bound procedures
!!#############################################################################
!! ========== Public Implementations ==========
!! --- Utilities ---
    Function ks_array_allocated(A) Result(res)
      Integer :: res
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Type(info_ks_point) :: this_ks
      res = 0
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        If(A%ks_points(A%iterator_value)%data%allocated()) res = res + 1
        this_ks = A%iterator_next()
      Enddo
    End Function ks_array_allocated

    Function ks_array_initialized(A) Result(res)
      Logical :: res
      Class(ks_array), Intent(InOut) :: A
      res = Allocated(A%ks_points)
    End Function ks_array_initialized

    Function ks_array_get_number_of_ks_points(A) Result(res)
    Use dash_utils, Only: dash_error
      Integer                     :: res
      Class(ks_array), Intent(In) :: A
      res = A%n_ks_points
      !! just a check
      If(Allocated(A%ks_points)) Then
        If(Size(A%ks_points).ne.res) Call dash_error(1,'ks_array_get_number_of_ks_points','size of ks_array not coherent')
      Else
        If(res.ne.0) Call dash_error(1,'ks_array_get_number_of_ks_points','size of ks_array not coherent')
      Endif
    End Function ks_array_get_number_of_ks_points

    Subroutine global_store_ks_array(A)
    !! Allocate "permanently" points in memory
      Class(ks_array), Intent(InOut) :: A
      ! Local variables
      Type(info_ks_point) :: this_ks
      !! Do not do anything if low memory option active
      !! Just a place holder for now
      If(memory_execution_mode.eq.MEMORY_EXECUTION_LOW) Return
      this_ks = A%iterator_init()
      Do While (this_ks%exist())
        If(.not.A%ks_points(A%iterator_value)%data%allocated()) Then
          Call A%ks_points(A%iterator_value)%data%store(1,this_ks%ks_size)
        Endif
        this_ks = A%iterator_next()
      Enddo
    End Subroutine global_store_ks_array

    Subroutine global_free_ks_array(A)
      Class(ks_array), Intent(InOut) :: A
      ! Local variables
      Type(info_ks_point) :: this_ks
      !! Do not do anything if low memory option active
      If(memory_execution_mode.eq.MEMORY_EXECUTION_LOW) Return
      this_ks = A%iterator_init()
      Do While (this_ks%exist())
        If(A%ks_points(A%iterator_value)%data%allocated()) Then
          !Call A%ks_points(A%iterator_value)%data%free(1,this_ks%ks_size)
          Call A%ks_points(A%iterator_value)%data%free()
        Else
          Write(*,*) "global_free_ks_array.... not allocated"
        Endif
        this_ks = A%iterator_next()
      Enddo
    End Subroutine global_free_ks_array

    Subroutine ks_array_set_indexing_from_array(A,raw_indeces)
    !! Set my_indeces in ks_array, if already resent some indeces add new ones after
    Use dash_utils
    Use interface_MPI
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
    End Subroutine ks_array_set_indexing_from_array

    Subroutine ks_array_print_distribution(A,weights,dist_aglo)
    Use interface_MPI!, Only: world
    Use backend_module, Only: bk_get_current_device
      Class(ks_array),               Intent(In)           :: A
      Real(double)   , Dimension(2), Intent(In), Optional :: weights
      Integer        ,               Intent(In), Optional :: dist_aglo
      ! Local variables
      Integer                                   :: my_rank
      Integer                                   :: i
      Integer                                   :: proc, max_proc
      Integer                                   :: ks, my_ks, max_n_ks
      Real(double)                              :: current_weight, max_weight, min_weight
      Integer     , Dimension(:)  , Allocatable :: my_ks_points
      Integer     , Dimension(:)  , Allocatable :: buffer
      Integer     , Dimension(:,:), Allocatable :: ks_data
      Integer     , Dimension(:)  , Allocatable :: buffer_devs
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
        Write(*,"(80('*'))")
        Write(*,"(A)", Advance='no') '       KS POINTS DISTRIBUTION'
        If(Present(dist_aglo)) Then
          Select Case(dist_aglo)
          Case(KS_DISTRIB_ALG_STATIC)
            Write(*,"(A)") ' STATIC'
          Case(KS_DISTRIB_ALG_COMPATIB)
            Write(*,"(A)") ' COMPATIBLE'
          Case(KS_DISTRIB_ALG_GREED)
            Write(*,"(A)") ' GREED'
          End Select
        Else
          Write(*,*)
        Endif
        max_weight = 0.0_double
        min_weight = Huge(min_weight)
        Do proc = 1, max_proc
          If(ks_data(1,proc).eq.INVALID) Cycle
          current_weight = 0.0_double
          i = 1
          Write(*,"('P',I5,': ')", Advance='no') proc - 1
          Do ks = 1, max_n_ks
            If(ks_data(ks,proc).eq.INVALID) Exit
            Associate( this_ks => A%all_points_info(ks_data(ks,proc)))
              Select Case(this_ks%ks_type)
              Case(POINT_IS_REAL)
                Write(*,"(A1,'(',I4,')')",Advance='no') 'R',this_ks%k_index
                If(Present(weights)) current_weight = current_weight + weights(1)
              Case(POINT_IS_COMPLEX)
                Write(*,"(A1,'(',I4,')')",Advance='no') 'C',this_ks%k_index
                If(Present(weights)) current_weight = current_weight + weights(2)
              End Select
              Select Case(this_ks%s_index)
              Case(SPIN_ALPHA)
                Write(*,"(A1,1X)", Advance='no') 'a'
              Case(SPIN_BETA)
                Write(*,"(A1,1X)", Advance='no') 'b'
              End Select
            End Associate
            If(i.ge.8) Then
              Write(*,*)
              Write(*,"(A8)",Advance='no') '      : '
              i = 0
            Endif
            i = i + 1
          Enddo
          Write(*,*)
          If(Present(weights)) Then
            Write(*,"(A15,F6.2)") '      : Weight ', current_weight
            max_weight = Max(max_weight,current_weight)
            min_weight = Min(min_weight,current_weight)
          Endif
          Write(*,"(A15,2X,I2)") '      : Device ', buffer_devs(proc)
        Enddo
        If(Present(weights)) Write(*,"(A,F10.2)") 'Maximum weight difference: ', max_weight - min_weight
        Deallocate(buffer)
        Deallocate(buffer_devs)
        Write(*,"(80('*'))")
      Endif
    End Subroutine ks_array_print_distribution

    Subroutine ks_array_create(A,source)
    !! Initialize the ks_array from a source
    !! Allocate requred space for the array of points ans sets all the info
    !! Do not allocate actual data for matrices
    Use dash_utils
      Class(ks_array), Intent(InOut) :: A
      Type(ks_array) , Intent(In   ) :: source
      ! Local variables
      Logical :: is_complex
      Integer :: nks
      Integer :: ks
      Integer :: i
      Integer :: ierr
      nks = source%n_ks_points      !! points in the mold array
      A%n_ks_points      = nks
      A%n_ks_points_type = source%n_ks_points_type
      A%all_points_info  = source%all_points_info
      !! If the array is already created exit
      !! This is necessary as we may call create also for array already created, inside SCF
      If( A%initialized() ) Return
      !! If the array does not contain any ks then retrun
      If(nks.eq.0) Return
      A%my_indeces = source%my_indeces !! Automatic Fortran reallocation on assignment
      !! Juast a Check
      If(Size(A%my_indeces).ne.nks) Call dash_error(0,'ks_array_create','wrong indeces count',Size(A%my_indeces))
      !! Allocate the required space for ks points
      Allocate(A%ks_points(nks), stat=ierr)
      If(ierr.ne.0) Call dash_error(0,'ks_array_create','allocation A%ks_points failed')
      Loop_on_owned_ks: Do ks = 1, nks
        !! Copy info from all to specific point_data
        A%ks_points(ks)%info = A%all_points_info(A%my_indeces(ks))
        !! Allocate the ks_point (ks_point_data%data) inside each ks_point_data
        is_complex = A%ks_points(ks)%info%ks_type.eq.POINT_IS_COMPLEX
        !! 1 is the number of irreps (NO SYMM: always 1...)
        Call A%ks_points(ks)%data%create(is_complex,1)
        A%ks_points(ks)%status = POINT_STATUS_INIT
        !! At this point the matrix is allocated according to the type of the k point
        !! (real or complex) but the actual data structure is not allocated.
      Enddo Loop_on_owned_ks
    End Subroutine ks_array_create

    Subroutine ks_array_create_vary(A,info_array)
    !! Generate a ks_array from a info_ks_point array
    !! The ks indeces of the array should be defined beforehand
    !! WARNING:: this procedure do not allocate the points vector!
    Use dash_utils
      Class(ks_array)    ,               Intent(InOut) :: A
      Type(info_ks_point), Dimension(:), Intent(In   ) :: info_array
      ! Local Variables
      Logical :: is_complex
      Integer :: nks
      Integer :: k, i
      If(.Not.Allocated(A%my_indeces)) Then
        A%all_points_info = info_array
        A%n_ks_points = 0
        A%n_ks_points_type(1) = 0
        A%n_ks_points_type(2) = 0
        Return
      Endif
      nks = Size(A%my_indeces)
      A%all_points_info = info_array
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
          A%n_ks_points = 0
          A%n_ks_points_type(1) = 0
          A%n_ks_points_type(2) = 0
        End Select
      Enddo
    End Subroutine ks_array_create_vary

!! --- LA Operations ---
    Subroutine ks_array_initialize_ortogonalization_matrices(A)
    !! Construction of ortogonalization matrices
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        Associate(iv => A%iterator_value)
          Call A%ks_points(iv)%data%cholesky()    !! get U from A = U^T * U
          Call A%ks_points(iv)%data%clean_lower() !! set lower triang. to 0
          Call A%ks_points(iv)%data%invert()      !! Compute inverted matrix
        End Associate
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_initialize_ortogonalization_matrices

    Subroutine ks_array_multiply(A,B,opb,C,opc)
    !! Multiplication of 2 ks_arrays
      Class(ks_array), Intent(InOut) :: A
      Type(ks_array) , Intent(InOut) :: B
      Integer        , Intent(In   ) :: opb
      Type(ks_array) , Intent(InOut) :: C
      Integer        , Intent(In   ) :: opc
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        Associate(iv => A%iterator_value)
          If(A%ks_points(iv)%data%allocated()) Then
            Call B%ks_points(iv)%data%dagger(opb)
            Call C%ks_points(iv)%data%dagger(opc)
            Call A%ks_points(iv)%data%multiply(B%ks_points(iv)%data,C%ks_points(iv)%data)
          Endif
        End Associate
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_multiply

    Subroutine ks_array_diag(A,evals)
      Class(ks_array)  , Intent(InOut) :: A
      Type(ks_array_1D), Intent(InOut) :: evals
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        Associate( Aks_dat => A%ks_points(A%iterator_value)%data, &
                   e_dativ => evals%data(A%iterator_value)         )
          Call Aks_dat%diag(e_dativ)
        End Associate
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_diag

    Subroutine ks_array_copy(A,B)
    !! Transfer arrays
      Class(ks_array), Intent(InOut) :: A
      Type(ks_array) , Intent(InOut) :: B
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        Associate(iv => A%iterator_value)
          !If(A%ks_points(iv)%data%allocated()) Then
            Call A%ks_points(iv)%data%copy(B%ks_points(iv)%data)
          !Endif
        End Associate
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_copy

!! --- Iterative Operations ---
    Subroutine ks_array_shift_diagonal(A,shift)
    !! Shift values on the main diagonal, only occupied ones
      Class(ks_array), Intent(InOut) :: A
      Real(double)   , Intent(In   ) :: shift
      ! Local Variables
      Type(info_ks_point) :: this_ks
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        If(this_ks%nocc.gt.0) Then
          Call A%ks_points(A%iterator_value)%data%shift_diagonal(shift,this_ks%nocc)
        Endif
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_shift_diagonal

!! --- General Operations ---
    Subroutine ks_array_set_nocc(A,jalph,n_k_points,n_spin_states)
    Use dash_utils, Only: dash_error
    !! Set the number of occupied states for the ks in this array
    !! This is somehow a hack, shold be improved...
      Class(ks_array) ,               Intent(InOut) :: A
      Integer         , Dimension(:), Intent(In   ) :: jalph
      Integer         ,               Intent(In   ) :: n_k_points
      Integer         ,               Intent(In   ) :: n_spin_states
      !! jalph is the jalpha from CRYSTAL, this contains the number of occupied
      !! orbitals for each k points.
      !! The order is K1S1 K1S2 K2S1 K2S2 K3S1 K3S2...
      ! Local Variables
      Type(info_ks_point) :: this_ks
      Integer, Dimension(1:n_k_points) :: k_occ_alpha
      Integer, Dimension(1:n_k_points) :: k_occ_beta
      Integer :: k
      Integer :: nks
      Integer :: rrr
      nks = n_k_points*n_spin_states
      If(nks.ne.Size(A%all_points_info)) Call dash_error(0,"set_nocc","Wrong array dimension")
      Select Case(n_spin_states)
      Case(CLOSE_SHELL)
        k_occ_alpha = jalph(1:nks)
      Case(OPEN_SHELL)
        k_occ_alpha = jalph(1:nks:2)
        k_occ_beta  = jalph(2:nks:2)
      End Select
      !! Set occupation level for all_points_info,
      !! based on k index and spin index
      Do k = 1, Size(A%all_points_info)
        Select Case(A%all_points_info(k)%get_spin())
        Case(SPIN_ALPHA)
          A%all_points_info(k)%nocc = k_occ_alpha(A%all_points_info(k)%get_k())
        Case(SPIN_BETA)
          A%all_points_info(k)%nocc = k_occ_beta(A%all_points_info(k)%get_k())
        End Select
      Enddo
      !! Set occupation level for the actual array of this process
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        Select Case(this_ks%get_spin())
        Case(SPIN_ALPHA)
          A%ks_points(A%iterator_value)%info%nocc = k_occ_alpha(this_ks%get_k())
        Case(SPIN_BETA)
          A%ks_points(A%iterator_value)%info%nocc = k_occ_beta(this_ks%get_k())
        End Select
        this_ks = A%iterator_next()
      Enddo
    End Subroutine ks_array_set_nocc

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
        Allocate(A1D%data(1:1))
        Allocate(A1D%info(1:1))
        A1D%info(1)%ks_type = POINT_NOT_EXIST
        A1D%status = POINT_STATUS_ALLOC
        Return
      Endif
      Allocate(A1D%data(1:nks))
      Allocate(A1D%info(1:nks))
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        Associate( iv => A%iterator_value )
          A1D%info(iv) = this_ks
          Call A1D%data(iv)%create(ldim)
        End Associate
        this_ks = A%iterator_next()
      Enddo
      A1D%status = POINT_STATUS_ALLOC
    End Subroutine ks_array_1D_create

    Subroutine ks_array_1D_destroy(A1D)
    !! Deallocate the space for 1D array
    !! Destroy all info inside
      Class(ks_array_1D), Intent(InOut) :: A1D
      ! Local Variables
      Integer :: ks
      Type(info_ks_point) :: this_ks
      Do ks = 1, Size(A1D%info)
        If(A1D%info(ks)%ks_type.ne.POINT_NOT_EXIST) Call A1D%data(ks)%destroy()
      Enddo
      Deallocate(A1D%data)
      Deallocate(A1D%info)
      A1D%status = POINT_STATUS_INIT
    End Subroutine ks_array_1D_destroy

    Subroutine ks_array_1D_get_raw(A1D,raw_data,ind)
      Class(ks_array_1D),                            Intent(InOut) :: A1D
      Real(double)      , Dimension(:), Allocatable, Intent(  Out) :: raw_data
      Integer           ,                            Intent(In   ) :: ind
      Call A1D%data(ind)%get_raw(raw_data)
    End Subroutine ks_array_1D_get_raw

    Subroutine ks_array_1D_set_raw(A1D,raw_data,ind)
      Class(ks_array_1D),               Intent(InOut) :: A1D
      Real(double)      , Dimension(:), Intent(In   ) :: raw_data
      Integer           ,               Intent(In   ) :: ind
      Call A1D%data(ind)%set_raw(raw_data)
    End Subroutine ks_array_1D_set_raw

    Subroutine ks_array_set_raw_real(A,raw_data,ind)
    ! Insert raw_data in ks point at position ind
      Class(ks_array),                 Intent(InOut) :: A
      Real(double)   , Dimension(:,:), Intent(In   ) :: raw_data
      Integer        ,                 Intent(In   ) :: ind

      Call A%ks_points(ind)%data%set_raw(raw_data)
    End Subroutine ks_array_set_raw_real
    Subroutine ks_array_set_raw_comp(A,raw_data,ind)
    ! Insert raw_data in ks point at position ind
      Class(ks_array),                 Intent(InOut) :: A
      Complex(double), Dimension(:,:), Intent(In   ) :: raw_data
      Integer        ,                 Intent(In   ) :: ind

      Call A%ks_points(ind)%data%set_raw(raw_data)
    End Subroutine ks_array_set_raw_comp

    Subroutine ks_array_get_raw_real(A,raw_data,ind)
    ! Extract raw_data in ks point at position ind
      Class(ks_array),                              Intent(InOut) :: A
      Real(double)   , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Integer        ,                              Intent(In   ) :: ind

      Call A%ks_points(ind)%data%get_raw(raw_data)
    End Subroutine ks_array_get_raw_real
    Subroutine ks_array_get_raw_vector(A,raw_vect,ind,k_types)
    ! Extract raw_data in ks point at position ind
      Class(ks_array),                            Intent(InOut) :: A
      Real(double)   , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
      Integer        ,                            Intent(In   ) :: ind
      Integer        ,                            Intent(In   ) :: k_types

      Call A%ks_points(ind)%data%get_raw(raw_vect,k_types)
    End Subroutine ks_array_get_raw_vector
    Subroutine ks_array_get_raw_comp(A,raw_data,ind)
    ! Extract raw_data in ks point at position ind
      Class(ks_array),                              Intent(InOut) :: A
      Complex(double), Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Integer        ,                              Intent(In   ) :: ind

      Call A%ks_points(ind)%data%get_raw(raw_data)
    End Subroutine ks_array_get_raw_comp

    Subroutine initialize_ks_array(A,nks)
    Use dash_utils
    !! Initialize the array based on how many points it should contain
      Class(ks_array), Intent(  Out) :: A
      Integer        , Intent(In   ) :: nks       !! Total number of ks points
      ! Local Variables
      Integer :: ierr
      If(Allocated(A%all_points_info)) Deallocate(A%all_points_info)
      Allocate(A%all_points_info(nks), stat=ierr)
      If(ierr.ne.0) Call dash_error(0,'initialize_ks_array','allocating all_points_info')
    End Subroutine initialize_ks_array

    Subroutine initialize_ks_array_from_source(A,source)
    Use dash_utils
    !! Initialize the array based on how many points it should contain
      Class(ks_array), Intent(  Out) :: A
      Type(ks_array) , Intent(In   ) :: source    !! mold array
      ! Local variables
      Integer :: n
      Integer :: ierr
      n = Size(source%all_points_info)
      If(Allocated(A%all_points_info)) Deallocate(A%all_points_info)
      Allocate(A%all_points_info(n), stat=ierr)
      If(ierr.ne.0) Call dash_error(0,'initialize_ks_array_from_source','allocating all_points_info')
    End Subroutine initialize_ks_array_from_source

    Function ks_array_iterator_next(A) Result(info)
      Class(ks_array)    , Intent(InOut) :: A
      Type(info_ks_point)                :: info

      !Allocate(info)
      A%iterator_value = A%iterator_value + 1
      If(A%iterator_value.le.Ubound(A%ks_points, Dim=1))Then
        info = A%ks_points(A%iterator_value)%info
      Else
        A%iterator_value = Ubound(A%ks_points, Dim=1) + 1
        info%ks_type = POINT_NOT_EXIST
      Endif
    End Function ks_array_iterator_next

    Function ks_array_iterator_prev(A) Result(info)
      Class(ks_array)    , Intent(InOut) :: A
      Type(info_ks_point)                :: info

      !Allocate(info)
      A%iterator_value = A%iterator_value - 1
      If(A%iterator_value.ge.Lbound(A%ks_points, Dim=1))Then
        info = A%ks_points(A%iterator_value)%info
      Else
        A%iterator_value = Lbound(A%ks_points, Dim=1) - 1
        info%ks_type = POINT_NOT_EXIST
      Endif
    End Function ks_array_iterator_prev

    Subroutine ks_array_iterator_reset(A)
      Class(ks_array), Intent(InOut) :: A

      A%iterator_value = INVALID
    End Subroutine ks_array_iterator_reset

    Function ks_array_iterator_init(A) Result(info)
      Class(ks_array)    , Intent(InOut) :: A
      Type(info_ks_point)                :: info

      A%iterator_value = 1
      If(.not.Allocated(A%ks_points))Then
        info%ks_type = POINT_NOT_EXIST
        A%iterator_value = INVALID
      Else
        info = A%ks_points(A%iterator_value)%info
      Endif
    End Function ks_array_iterator_init

    Subroutine local_store_ks_array(A)
      Class(ks_array), Intent(InOut) :: A
      !! NOT IMPLEMENTED
    End Subroutine local_store_ks_array

    Subroutine local_free_ks_array(A)
      Class(ks_array), Intent(InOut) :: A
      !! NOT IMPLEMENTED
    End Subroutine local_free_ks_array

    Function ks_array_compute_tester(A) Result (tester)
    Use base_numbers
    !! Compute the tester value, used to check SCF convergence
      Real(double)                   :: tester
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Real(double)        :: tmp_tester
      Type(info_ks_point) :: this_ks

      tester = 0.0_double
      this_ks = A%iterator_init()
      Do While(this_ks%exist())
        If(this_ks%nocc.gt.0) Then
          tmp_tester = A%ks_points(A%iterator_value)%data%compute_tester(this_ks%nocc)
          tester = Max(tmp_tester,tester)
        Endif
        this_ks = A%iterator_next()
      Enddo
    End Function ks_array_compute_tester

    Logical Function info_ks_point_exist(I)
    Use dash_utils, Only: dash_error
      Class(info_ks_point), Intent(InOut) :: I
      Selectcase(I%ks_type)
      Case(POINT_NOT_EXIST)
        info_ks_point_exist = .False.
      Case(POINT_IS_REAL)
        info_ks_point_exist = .True.
      Case(POINT_IS_COMPLEX)
        info_ks_point_exist = .True.
      Case Default
        info_ks_point_exist = .False.
        Call dash_error(0,'info_ks_point_exist','invalid k type')
      Endselect
    End Function

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

    Integer Function info_ks_point_get_nocc(I)
    Class(info_ks_point), Intent(In) :: I
    info_ks_point_get_nocc = I%nocc
    End Function info_ks_point_get_nocc

    Subroutine ks_array_hack(A, res)
    ! Construction of ortogonalization matrices
    Use base_numbers
      Type(ks_array) , Intent(  Out) :: res
      Class(ks_array), Intent(InOut) :: A
      ! Local Variables
      Type(info_ks_point) :: this_ks
      Real(double), Dimension(:,:), Allocatable :: buffer
      Integer :: i
      Call res%create(A)               !! Allocating space for the result
      Call Move_alloc(res%ks_points, A%ks_points)
    End Subroutine ks_array_hack

    Subroutine ks_array_move(lhs,rhs)
      Class(ks_array), Intent(  Out) :: lhs
      Class(ks_array), Intent(In   ) :: rhs
      ! Local Variables
      Type(info_ks_point) :: this_ks
      Call lhs%create(rhs)
      lhs%ks_points = rhs%ks_points
    End Subroutine ks_array_move

    Subroutine Print_ks_array_all_points_info(A, out, name)
      Class(ks_array),           Intent(In) :: A
      Integer            , Optional, Intent(In) :: out
      Character(Len=*)   , Optional, Intent(In) :: name
      ! Local Variables
      Integer :: iuni
      Integer :: k
      If(Present(out))Then
        iuni = out
      Else
        iuni = 6
      Endif
      If(Present(name))Then
        Write(iuni,"(1X,A,' contains ',I0,' ks points:')") name, A%n_ks_points
      Else
        Write(iuni,"(1X,'ks_array contains ',I0,' ks points:')") A%n_ks_points
      Endif
    End Subroutine Print_ks_array_all_points_info

    Subroutine increment_ks_points_counter(A,point_type)
    ! Update the counter for the number and type of points
      Class(ks_array), Intent(InOut) :: A
      Integer        , Intent(In   ) :: point_type
      Select case(point_type)
      Case(POINT_IS_COMPLEX)
        A%n_ks_points_type(2) = A%n_ks_points_type(2) + 1
        A%n_ks_points         = A%n_ks_points         + 1
      Case(POINT_IS_REAL   )
        A%n_ks_points_type(1) = A%n_ks_points_type(1) + 1
        A%n_ks_points         = A%n_ks_points         + 1
      Case Default
        Write(*,*) 'ERROR:: Invalid point type to assign'
        Stop
      End Select
    End Subroutine increment_ks_points_counter

End Module module_arrays
