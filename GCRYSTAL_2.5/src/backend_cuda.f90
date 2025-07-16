Module backend_module
Use base_numbers
Use cudafor
Use cublas_v2
Use cusolverdn

  Implicit None
  Private

  Logical, Private :: do_HOST_lapack = .False.
  Logical, Private :: do_DEV_cuda    = .True.

  Integer, Parameter, Private :: MAX_THREADS    = 1024
  Integer, Parameter, Private :: MAX_THREADS_2D = 32

  !! --- VARIABLES ---
  Integer               , Private :: assigned_dev = INVALID   !! only processes with gpus performs CUDA operations
  Type(cublashandle)    , Private :: blas_handle              !! CUDA stuff
  Type(cusolverdnhandle), Private :: solver_handle
  Type(cusolverdnparams), Private :: solver_params
  Integer               , Private :: init_blas_handle_status   = 0
  Integer               , Private :: init_solver_handle_status = 0
  Integer               , Private :: init_solver_params_status = 0
  Integer               , Public  :: bk_magic_factor           = 1

  !! --- DERIVED TYPES ---
  Type, Public :: bk_real_matrix
    Real(fp64), Dimension(:,:), Allocatable         :: data
    Real(fp64), Dimension(:,:), Device, Allocatable :: data_d
    Integer                                         :: n = INVALID
    Integer                                         :: m = INVALID
  End Type bk_real_matrix
  Type, Public :: bk_comp_matrix
    Complex(fp64), Dimension(:,:), Allocatable         :: data
    Complex(fp64), Dimension(:,:), Device, Allocatable :: data_d
    Integer                                            :: n = INVALID
    Integer                                            :: m = INVALID
  End Type bk_comp_matrix
  Type, Public :: bk_vector_1D
    Real(fp64), Dimension(:), Allocatable         :: data
    Real(fp64), Dimension(:), Allocatable, Device :: data_d
    Integer                                       :: n = INVALID
  End Type bk_vector_1D

  Type, Public :: ft_data
    Integer, Dimension(3)    ,              Private          :: k_grid
    Integer, Dimension(3)    ,              Private          :: k_grid_mult
    Integer,                                Private          :: k_grid_lcm
    Integer,                                Private          :: n_couples
    Integer, Dimension(:)    , Allocatable, Private, Managed :: shell1_size
    Integer, Dimension(:)    , Allocatable, Private, Managed :: shell1_start
    Integer, Dimension(:)    , Allocatable, Private, Managed :: shell2_size
    Integer, Dimension(:)    , Allocatable, Private, Managed :: shell2_start
    Integer, Dimension(:)    , Allocatable, Private, Managed :: pg_shell_start
    Integer, Dimension(:)    , Allocatable, Private, Managed :: repetitions
    Integer, Dimension(:,:,:), Allocatable, Private, Managed :: g_vecs
    Integer,                                Private          :: status = STATUS_INVALID
    Contains
    Procedure, Public  :: set   => ft_data_set
    Procedure, Public  :: free  => ft_data_free
    Procedure, Public  :: print => ft_data_print
    Procedure, Private :: size  => ft_data_get_alloc_mem
  End type
  !! Data for constructing direct-space vector with CRYSTAL memory layout
  Type(ft_data), Public :: ft_map

  !! --- PREPROCESSING MACROS ---
#define cuda_safe_Call( err ) call __CUDA_safe_call( err, 'backend_cuda.f90', __LINE__ )
#define cuda_error( msg ) call __CUDA_error_msg( msg, 'backend_cuda.f90', __LINE__ )
#define cuda_check( err ) call __CUDA_error_check( err, 'backend_cuda.f90', __LINE__ )
#define culib_check( err ) call __CUDA_lib_error_check( err, 'backend_cuda.f90', __LINE__ )
#define cuda_get_error_check call __CUDA_error_check_var( 'backend_cuda.f90', __LINE__ )
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef DO_MEM_PRINT
  #define print_device_memory( name ) call bk_mem_print(name)
#else
  #define print_device_memory( name )
#endif
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - -

  !! --- PROCEDURES ---
  !! --- Utility ---
  Public :: bk_get_backend
  Public :: bk_get_visible_devices
  Public :: bk_get_current_device
  Public :: bk_set_current_device
  Public :: bk_initialize_gpus
  Public :: bk_mem_print
  Public :: bk_get_memory_device
  !! --- Memory Managment ---
  Public :: bk_real_alloc  , bk_comp_alloc
  Public :: bk_real_dealloc, bk_comp_dealloc
  Public :: bk_alloc_1d
  Public :: bk_dealloc_1d
  !! --- Data Movement ---
  Public :: bk_set_raw_real, bk_set_raw_comp
  Public :: bk_get_raw_real, bk_get_raw_comp
  Public :: bk_set_raw_1d
  Public :: bk_get_raw_1d
  !! --- LA Operations ---
  Public :: bk_real_cholesky       , bk_comp_cholesky
  Public :: bk_real_set_lower_to   , bk_comp_set_lower_to
  Public :: bk_real_invert         , bk_comp_invert
  Public :: bk_real_mult           , bk_comp_mult
  Public :: bk_real_sym_diag       , bk_comp_sym_diag
  Public :: bk_real_matrix_copy    , bk_comp_matrix_copy

  Public :: bk_sort_eval_evecs
  Interface bk_sort_eval_evecs
    Module Procedure :: bk_real_sort_eval_evecs
    Module Procedure :: bk_comp_sort_eval_evecs
  End Interface bk_sort_eval_evecs


  !! --- Density matrix ---
  Public :: bk_real_matrix_compute_density, bk_comp_matrix_compute_density
  !! --- Matrix Related ---
  Public :: bk_real_allocated         , bk_comp_allocated
  Public :: bk_real_shift_diag        , bk_comp_shift_diag
  Public :: bk_compute_real_tester_sub, bk_compute_complex_tester_sub
  Public :: bk_shift_1d
  !! --- Fermi ---
  !Public :: bk_setup_fermi_data
  !Public :: bk_get_bands_range_1d
  !Public :: bk_get_fermi_data
  !Public :: bk_get_fermi_ranges
  !! --- Host Memory Storage ---
  Public :: bk_host_real_alloc      , bk_host_comp_alloc
  Public :: bk_host_real_dealloc    , bk_host_comp_dealloc
  Public :: bk_host_real_matrix_copy, bk_host_comp_matrix_copy
  Public :: bk_host_set_raw_real    , bk_host_set_raw_comp
  Public :: bk_host_get_raw_real    , bk_host_get_raw_comp
  Public :: bk_real_host_to_dev     , bk_comp_host_to_dev
  Public :: bk_real_dev_to_host     , bk_comp_dev_to_host

  Contains

!! --- Utility ---
    Integer Function bk_get_backend()
    !! return backend identifier index
      bk_get_backend = 2
    End Function bk_get_backend

    Integer Function bk_get_visible_devices(limit)
    !! - DEVICE QUERY -
    !! return number of available decices for this process
    !! never higher than limit
      Integer, Intent(In), Optional :: limit
      ! Local Variables
      Integer :: n_dev
      cuda_safe_Call( cudaGetDeviceCount(n_dev) )
      If(Present(limit)) Then
        If(n_dev.gt.limit) n_dev = limit
      Endif
      bk_get_visible_devices = n_dev
    End Function bk_get_visible_devices

    Integer Function bk_get_current_device()
    !! - DEVICE QUERY -
    !! return current device associated with this process
      bk_get_current_device = INVALID
      If(query_check()) cuda_safe_Call( cudaGetDevice(bk_get_current_device) )
    End Function bk_get_current_device

    Subroutine bk_set_current_device(index,limit)
    !! - DEVICE QUERY -
    !! set gpu index for thsi process, if index too high set gpu 0
      Integer, Intent(In)           :: index
      Integer, Intent(In), Optional :: limit
      ! Local Variables
      Integer :: ierr, avail_devs, idx
      cuda_safe_Call( cudaGetDeviceCount(avail_devs) )
      If(Present(limit)) Then
        If(avail_devs.gt.limit) avail_devs = limit
      Endif
      If(index.gt.avail_devs) Then
        idx = 0
        cuda_error("required device index too high, device 0 used instead")
      Else
        idx = index
      Endif
      cuda_safe_Call( cudaSetDevice(idx) )
      assigned_dev = idx          !! now this porcess has a device
    End Subroutine bk_set_current_device

    Subroutine bk_initialize_gpus()
    !! - DEVICE OPERATION -
    !! Initialize cuSOLVER and cuBLAS handles for this process in current assigned device
      Integer :: ierr
      If(query_check()) Then
        If(init_blas_handle_status.eq.0) Then
          cuda_safe_Call( cublascreate(blas_handle) )
          init_blas_handle_status = 1
        Endif
        If(init_solver_handle_status.eq.0) Then
          cuda_safe_Call( cusolverdncreate(solver_handle) )
          init_solver_handle_status = 1
        Endif
        If(init_solver_params_status.eq.0) Then
          cuda_safe_Call( cusolverdncreateparams(solver_params) )
          init_solver_params_status = 1
        Endif
      Else
        cuda_error("initializing gpu handles when devices are not assigned")
      Endif
    End Subroutine bk_initialize_gpus

    Subroutine bk_mem_print(name)
    !! print in output memory usage of GPU currently associated with this process
    !! ONLY USED FOR DEBUG PURPOSES
      Character(Len=*), Intent(In) :: name
      ! Local Variables
      Integer                       :: ierr, dev
      Integer(kind=cuda_count_kind) :: free, total
      Real(fp64)                  :: free_GB, total_GB
      cuda_safe_Call( cudaGetDevice(dev) )
      cuda_safe_Call( cudaMemGetInfo(free,total) )
      free_GB  = Real(free,Kind=fp64)/(1024.0_fp64**3)
      total_GB = Real(total,Kind=fp64)/(1024.0_fp64**3)
      Write(*,"(A,T26,A,I1,A,F7.3,A,F7.3,A)") name," =======CUDA=====> DEVICE ",dev," MEMORY:",total_GB-free_GB," / ",total_GB," GB"
    End Subroutine bk_mem_print

    Subroutine bk_get_memory_device(used_mem,total_mem)
    !! - DEVICE QUERY -
    !! return used memory and total memory (MB) of the GPU currently
    !! associated with this process
      Real(fp64), Intent(Out) :: used_mem
      Real(fp64), Intent(Out) :: total_mem
      ! Local Variables
      Integer(Kind=cuda_count_kind) :: free_byte, total_byte
      free_byte = Int(0,Kind=cuda_count_kind)
      total_byte = Int(0,Kind=cuda_count_kind)
      If(query_check()) cuda_safe_Call( cudaMemGetInfo(free_byte,total_byte) )
      total_mem = Real(total_byte,Kind=fp64)/(1024.0_fp64**2)
      used_mem = total_mem - Real(free_byte,Kind=fp64)/(1024.0_fp64**2)
    End Subroutine bk_get_memory_device

!! --- Memory Managment ---
    Subroutine bk_real_alloc(a,n,m)
      Type(bk_real_matrix), Intent(InOut) :: a
      Integer             , Intent(In   ) :: n
      Integer             , Intent(In   ) :: m
      ! Local Variables
      Integer :: ierr
      a%n = n
      a%m = m
      Allocate(a%data_d(n,m), stat=ierr)
      cuda_check(ierr)
      print_device_memory('bk_real_alloc')
    End Subroutine bk_real_alloc

    Subroutine bk_comp_alloc(a,n,m)
      Type(bk_comp_matrix), Intent(InOut) :: a
      Integer             , Intent(In   ) :: n
      Integer             , Intent(In   ) :: m
      ! Local Variables
      Integer :: ierr
      a%n = n
      a%m = m
      Allocate(a%data_d(n,m), stat=ierr)
      cuda_check(ierr)
      print_device_memory('bk_comp_alloc')
    End Subroutine bk_comp_alloc

    Subroutine bk_real_dealloc(a)
      Type(bk_real_matrix), Intent(InOut) :: a
      ! Local Variables
      Integer :: ierr
      If(.not.Allocated(a%data_d)) cuda_error("deallocating not-allocated-variable")
      Deallocate(a%data_d, stat=ierr)
      cuda_check(ierr)
      a%n = INVALID
      a%m = INVALID
      print_device_memory('bk_real_dealloc')
    End Subroutine bk_real_dealloc

    Subroutine bk_comp_dealloc(a)
      Type(bk_comp_matrix), Intent(InOut) :: a
      ! Local Variables
      Integer :: ierr
      If(.not.Allocated(a%data_d)) cuda_error("deallocating not-allocated-variable")
      Deallocate(a%data_d, stat=ierr)
      cuda_check(ierr)
      a%n = INVALID
      a%m = INVALID
      print_device_memory('bk_comp_dealloc')
    End Subroutine bk_comp_dealloc

    Subroutine bk_alloc_1d(vec,n)
      Type(bk_vector_1D), Intent(InOut) :: vec
      Integer           , Intent(In   ) :: n
      ! Local Variables
      Integer :: ierr
      If(Allocated(vec%data_d)) cuda_error("allocating alredy allocated variable")
      vec%n = n
      Allocate(vec%data_d(n), stat=ierr)
      cuda_check(ierr)
      print_device_memory('bk_alloc_1d')
    End Subroutine bk_alloc_1d

    Subroutine bk_dealloc_1d(vec)
      Type(bk_vector_1D), Intent(InOut) :: vec
      ! Local Variables
      Integer :: ierr
      If(.not.Allocated(vec%data_d)) cuda_error("allocating alredy allocated variable")
      Deallocate(vec%data_d, stat=ierr)
      cuda_check(ierr)
      vec%n = INVALID
      print_device_memory('bk_dealloc_1d')
    End Subroutine bk_dealloc_1d

!! --- Data Movement ---
    Subroutine bk_set_raw_real(a,raw_data)
      Type(bk_real_matrix),                 Intent(InOut) :: a
      Real(fp64)          , Dimension(:,:), Intent(In   ) :: raw_data
      If(Size(raw_data,Dim=1).ne.a%n) cuda_error("wrong domension")
      If(Size(raw_data,Dim=2).ne.a%m) cuda_error("wrong dimension")
      a%data_d = raw_data
      cuda_get_error_check
    End Subroutine bk_set_raw_real

    Subroutine bk_set_raw_comp(a,raw_data)
      Type(bk_comp_matrix),                 Intent(InOut) :: a
      Complex(fp64)       , Dimension(:,:), Intent(In   ) :: raw_data
      If(Size(raw_data,Dim=1).ne.a%n) cuda_error("wrong domension")
      If(Size(raw_data,Dim=2).ne.a%m) cuda_error("wrong dimension")
      a%data_d = raw_data
      cuda_get_error_check
    End Subroutine bk_set_raw_comp

    Subroutine bk_get_raw_real(a,raw_data)
      Type(bk_real_matrix),                              Intent(In   ) :: a
      Real(fp64)          , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      ! Local Variables
      Integer :: ldn, ldm
      ldn = a%n
      ldm = a%m
      Allocate(raw_data(ldn,ldm))
      raw_data = a%data_d
      cuda_get_error_check
    End Subroutine bk_get_raw_real

    Subroutine bk_get_raw_comp(a,raw_data)
      Type(bk_comp_matrix),                              Intent(In   ) :: a
      Complex(fp64)       , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      ! Local Variables
      Integer :: ldn, ldm
      ldn = a%n
      ldm = a%m
      Allocate(raw_data(ldn,ldm))
      raw_data = a%data_d
      cuda_get_error_check
    End Subroutine bk_get_raw_comp

    Subroutine bk_set_raw_1d(vec,raw)
      Type(bk_vector_1D),               Intent(InOut) :: vec
      Real(fp64)        , Dimension(:), Intent(In   ) :: raw
      If(Size(raw,Dim=1).ne.vec%n) cuda_error("wrong domension")
      vec%data_d = raw
      cuda_get_error_check
    End Subroutine bk_set_raw_1d

    Subroutine bk_get_raw_1d(vec,raw)
      Type(bk_vector_1D),                            Intent(InOut) :: vec
      Real(fp64)        , Dimension(:), Allocatable, Intent(  Out) :: raw
      ! Local Variables
      Integer :: ldn
      ldn = vec%n
      Allocate(raw(ldn))
      raw = vec%data_d
      cuda_get_error_check
    End Subroutine bk_get_raw_1d

!! --- LA Operations ---
    Subroutine bk_real_cholesky(matrix)
      Type(bk_real_matrix), Intent(InOut) :: matrix
      ! Local Variables
      Integer :: n
      Integer :: ierr
      Integer :: lwork
      ! Device Variables
      Real(fp64), Dimension(:), Allocatable, Device :: work_d
      Integer   ,                            Device :: ierr_d
      n = matrix%n
      call __CUDA_nvtx_profiler_start("real_cholesky_full",1)
      call __CUDA_nvtx_profiler_start("real_cholesky_buffer",2)
      cuda_safe_Call( cusolverdndpotrf_buffersize(solver_handle,CUBLAS_FILL_MODE_UPPER,n,matrix%data_d,n,lwork) )
      call __CUDA_nvtx_profiler_end
      Allocate(work_d(lwork),stat = ierr)
      cuda_check(ierr)
      call __CUDA_nvtx_profiler_start("real_cholesky_compute",2)
      cuda_safe_Call( cusolverdndpotrf(solver_handle,CUBLAS_FILL_MODE_UPPER,n,matrix%data_d,n,work_d,lwork,ierr_d) )
      call __CUDA_nvtx_profiler_end
      ierr = ierr_d
      culib_check(ierr)
      call __CUDA_nvtx_profiler_end
      print_device_memory('bk_real_cholesky')
    End Subroutine bk_real_cholesky

    Subroutine bk_comp_cholesky(matrix)
      Type(bk_comp_matrix), Intent(InOut) :: matrix
      ! Local Variables
      Integer :: n
      Integer :: ierr
      Integer :: lwork
      ! Device Variables
      Complex(fp64), Dimension(:), Allocatable, Device :: work_d
      Integer      ,                            Device :: ierr_d
      n = matrix%n
      cuda_safe_Call( cusolverdnzpotrf_buffersize(solver_handle,CUBLAS_FILL_MODE_UPPER,n,matrix%data_d,n,lwork) )
      Allocate(work_d(lwork),stat = ierr)
      cuda_check(ierr)
      cuda_safe_Call( cusolverdnzpotrf(solver_handle,CUBLAS_FILL_MODE_UPPER,n,matrix%data_d,n,work_d,lwork,ierr_d) )
      ierr = ierr_d
      culib_check(ierr)
      print_device_memory('bk_comp_cholesky')
    End Subroutine bk_comp_cholesky

    Subroutine bk_real_set_lower_to(matrix,val)
      Type(bk_real_matrix), Intent(InOut) :: matrix
      Real(fp64)          , Intent(In   ) :: val
      ! Local Variables
      Type(dim3) :: grid_size, block_size
      block_size = dim3(32,32,1)
      grid_size%x = (matrix%n + block_size%x-1) / block_size%x
      grid_size%y = (matrix%m + block_size%y-1) / block_size%y
      grid_size%z = 1
      call __CUDA_nvtx_profiler_start("real_set_lower_v2",1)
      Call real_set_lower_kernel_v2<<<grid_size,block_size,0,0>>>(matrix%data_d,val,matrix%n,matrix%m)
      cuda_get_error_check
      call __CUDA_nvtx_profiler_end
      print_device_memory('bk_real_set_lower_to')
    End Subroutine bk_real_set_lower_to

    Subroutine bk_comp_set_lower_to(matrix,val)
      Type(bk_comp_matrix), Intent(InOut) :: matrix
      Complex(fp64)       , Intent(In   ) :: val
      ! Local Variables
      !Integer :: n, i, j
      !Complex(fp64), Dimension(:,:), Device, Allocatable :: temp
      Type(dim3) :: grid_size, block_size
      block_size = dim3(32,32,1)
      grid_size%x = (matrix%n + block_size%x-1) / block_size%x
      grid_size%y = (matrix%m + block_size%y-1) / block_size%y
      grid_size%z = 1
      call __CUDA_nvtx_profiler_start("comp_set_lower_v2",1)
      Call comp_set_lower_kernel<<<grid_size,block_size,0,0>>>(matrix%data_d,val,matrix%n,matrix%m)
      cuda_get_error_check
      call __CUDA_nvtx_profiler_end
      print_device_memory('bk_comp_set_lower_to')
    End Subroutine bk_comp_set_lower_to

    Subroutine bk_real_invert(matrix)
      Type(bk_real_matrix), Intent(InOut) :: matrix
      ! Local Variables
      Integer                                   :: n
      Integer                                   :: ierr
      Integer(fp64)                             :: dlwork, hlwork
      Integer(int64), Parameter                 :: ONE = 1_int64
      Integer(int8) , Dimension(:), Allocatable :: hwork
      ! Device Variables
      Integer(int8), Dimension(:), Allocatable, Device :: dwork_d
      Integer      ,                            Device :: istat_d
      n = matrix%n
      call __CUDA_nvtx_profiler_start("real_invert_full",1)
      call __CUDA_nvtx_profiler_start("real_invert_buffer",2)
      cuda_safe_Call( cusolverdnxtrtri_buffersize(solver_handle,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,n,CUDADATATYPE(CUDA_R_64F),matrix%data_d,n,dlwork,hlwork) )
      call __CUDA_nvtx_profiler_end
      dlwork = Max(ONE,dlwork) * bk_magic_factor
      Allocate(dwork_d(dlwork),stat=ierr)
      cuda_check(ierr)
      hlwork = Max(ONE,hlwork) * bk_magic_factor
      Allocate(hwork(hlwork),stat=ierr)
      call __CUDA_nvtx_profiler_start("real_invert_compute",2)
      cuda_safe_Call( cusolverdnxtrtri(solver_handle,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,n,CUDADATATYPE(CUDA_R_64F),matrix%data_d,n,dwork_d,dlwork,hwork,hlwork,istat_d) )
      call __CUDA_nvtx_profiler_end
      ierr = istat_d
      cuda_check(ierr)
      call __CUDA_nvtx_profiler_end
      print_device_memory('bk_real_invert')
    End Subroutine bk_real_invert

    Subroutine bk_comp_invert(matrix)
      Type(bk_comp_matrix), Intent(InOut) :: matrix
      ! Local Variables
      Integer                                   :: n
      Integer                                   :: ierr
      Integer(fp64)                             :: dlwork, hlwork
      Integer(int64), Parameter                 :: ONE = 1_int64
      Integer(int8) , Dimension(:), Allocatable :: hwork
      ! Device Variables
      Integer(int8), Dimension(:), Allocatable, Device :: dwork_d
      Integer      ,                            Device :: istat_d
      n = matrix%n
      cuda_safe_Call( cusolverdnxtrtri_buffersize(solver_handle,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,n,cudadatatype(CUDA_C_64F),matrix%data_d,n,dlwork,hlwork) )
      dlwork = Max(ONE,dlwork) * bk_magic_factor
      Allocate(dwork_d(dlwork),stat=ierr)
      cuda_check(ierr)
      hlwork = Max(ONE,hlwork) * bk_magic_factor
      Allocate(hwork(hlwork),stat=ierr)
      cuda_safe_Call( cusolverdnxtrtri(solver_handle,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,n,cudadatatype(CUDA_C_64F),matrix%data_d,n,dwork_d,dlwork,hwork,hlwork,istat_d) )
      ierr = istat_d
      culib_check(ierr)
      print_device_memory('bk_comp_invert')
    End Subroutine bk_comp_invert

    Subroutine bk_real_mult(res,a,lopa,b,lopb)
    !! compurte res = alpha * op(a) * op(b) + beta * res
      Type(bk_real_matrix), Intent(InOut) :: res
      Type(bk_real_matrix), Intent(In   ) :: a
      Logical             , Intent(In   ) :: lopa
      Type(bk_real_matrix), Intent(In   ) :: b
      Logical             , Intent(In   ) :: lopb
      ! Local Variables
      Integer               :: ndim
      Character(Len=1)      :: opa
      Character(Len=1)      :: opb
      Integer               :: op_a
      Integer               :: op_b
      Real(fp64), Parameter :: ALPHA = 1.0_fp64
      Real(fp64), Parameter :: BETA  = 0.0_fp64
      Integer               :: ierr
      ndim = Max(a%n,b%n,res%n)
      op_a = Merge(CUBLAS_OP_T,CUBLAS_OP_N,lopa)
      op_b = Merge(CUBLAS_OP_T,CUBLAS_OP_N,lopb)
      cuda_safe_Call( cublasdgemm(blas_handle,op_a,op_b,ndim,ndim,ndim,ALPHA,a%data_d,ndim,b%data_d,ndim,BETA,res%data_d,ndim) )
      print_device_memory('bk_real_mult')
    End Subroutine bk_real_mult

    Subroutine bk_comp_mult(res,a,lopa,b,lopb)
    !! compurte res = alpha * op(a) * op(b)   +   beta * res
      Type(bk_comp_matrix), Intent(InOut) :: res
      Type(bk_comp_matrix), Intent(In   ) :: a
      Logical             , Intent(In   ) :: lopa
      Type(bk_comp_matrix), Intent(In   ) :: b
      Logical             , Intent(In   ) :: lopb
      ! Local Variables
      Integer                  :: ndim
      Character(Len=1)         :: opa
      Character(Len=1)         :: opb
      Integer                  :: op_a
      Integer                  :: op_b
      Complex(fp64), Parameter :: ALPHA = (1.0_fp64,0.0_fp64)
      Complex(fp64), Parameter :: BETA  = (0.0_fp64,0.0_fp64)
      Integer                  :: ierr
      ndim = Max(a%n,b%n,res%n)
      op_a = Merge(CUBLAS_OP_C,CUBLAS_OP_N,lopa)
      op_b = Merge(CUBLAS_OP_C,CUBLAS_OP_N,lopb)
      cuda_safe_Call( cublaszgemm(blas_handle,op_a,op_b,ndim,ndim,ndim,ALPHA,a%data_d,ndim,b%data_d,ndim,BETA,res%data_d,ndim) )
      print_device_memory('bk_comp_mult')
    End Subroutine bk_comp_mult

    Subroutine bk_real_sym_diag(matrix,evals)
    !! Compute all eigenvalues and eigenvectors of matrix
    !! matrix is lower triangular
      Type(bk_real_matrix), Intent(InOut) :: matrix
      Type(bk_vector_1D)  , Intent(InOut) :: evals
      ! Local Variables
      Integer                                   :: ierr
      Integer                                   :: ldm
      Integer(int64)                            :: hlwork, dlwork
      Integer(int8) , Dimension(:), Allocatable :: hwork
      ! Device Variables
      Integer(int8), Dimension(:), Allocatable, Device :: dwork_d
      Integer      ,                            Device :: istat_d
      ldm = matrix%n
      hlwork = -1_int64
      dlwork = -1_int64
      cuda_safe_Call( cusolverdnsetadvoptions(solver_params,0,CUSOLVER_ALG_0) )
      cuda_safe_Call( cusolverdnxsyevd_buffersize(solver_handle,solver_params,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,ldm, cudadatatype(CUDA_R_64F),matrix%data_d,ldm, cudadatatype(CUDA_R_64F),evals%data_d,cudadatatype(CUDA_R_64F),dlwork,hlwork) )
      Allocate(dwork_d(dlwork),stat=ierr)
      cuda_check(ierr)
      If(hlwork.le.0) hlwork = 1_int64
      Allocate(hwork(hlwork),stat=ierr)
      cuda_safe_Call( cusolverdnxsyevd(solver_handle,solver_params,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,ldm, cudadatatype(CUDA_R_64F),matrix%data_d,ldm, cudadatatype(CUDA_R_64F),evals%data_d,cudadatatype(CUDA_R_64F),dwork_d,dlwork,hwork,hlwork,istat_d) )
      ierr = istat_d
      culib_check(ierr)
      Deallocate(dwork_d)
      Deallocate(hwork)
      print_device_memory('bk_real_sym_diag')
    End Subroutine bk_real_sym_diag

    Subroutine bk_comp_sym_diag(matrix,evals)
    !! Compute all eigenvalues and eigenvectors of matrix
    !! matrix is lower triangular
      Type(bk_comp_matrix), Intent(InOut) :: matrix
      Type(bk_vector_1D)  , Intent(InOut) :: evals
      ! Local Variables
      Integer                                   :: ierr
      Integer                                   :: ldm
      Integer(int64)                            :: hlwork, dlwork
      Integer(int8) , Dimension(:), Allocatable :: hwork
      ! Device Variables
      Integer(int8), Dimension(:), Allocatable, Device :: dwork_d
      Integer      ,                            Device :: istat_d
      ldm = matrix%n
      hlwork = -1_int64
      dlwork = -1_int64
      cuda_safe_Call( cusolverdnsetadvoptions(solver_params,0,CUSOLVER_ALG_0) )
      cuda_safe_Call( cusolverdnxsyevd_buffersize(solver_handle,solver_params,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,ldm, cudadatatype(CUDA_C_64F),matrix%data_d,ldm,cudadatatype(CUDA_R_64F),evals%data_d,cudadatatype(CUDA_C_64F),dlwork,hlwork) )
      Allocate(dwork_d(dlwork),stat=ierr)
      cuda_check(ierr)
      If(hlwork.le.0) hlwork = 1_int64
      Allocate(hwork(hlwork),stat=ierr)
      cuda_safe_Call( cusolverdnxsyevd(solver_handle,solver_params,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,ldm,cudadatatype(CUDA_C_64F),matrix%data_d,ldm, cudadatatype(CUDA_R_64F),evals%data_d,cudadatatype(CUDA_C_64F),dwork_d,dlwork,hwork,hlwork,istat_d) )
      ierr = istat_d
      culib_check(ierr)
      Deallocate(dwork_d)
      Deallocate(hwork)
      print_device_memory('bk_comp_sym_diag')
    End Subroutine bk_comp_sym_diag

    Subroutine bk_real_matrix_copy(a,b)
      Type(bk_real_matrix), Intent(InOut) :: a
      Type(bk_real_matrix), Intent(InOut) :: b
      ! Local Variables
      Integer :: ierr
      If(.not.Allocated(a%data_d)) Then
        Allocate(a%data_d(b%n,b%m),stat=ierr)
        cuda_check(ierr)
      Endif
      a%data_d = b%data_d
      Deallocate(b%data_d,stat=ierr)
      cuda_check(ierr)
      b%n = INVALID
      b%m = INVALID
      print_device_memory('bk_real_matrix_copy')
    End Subroutine bk_real_matrix_copy

    Subroutine bk_comp_matrix_copy(a,b)
      Type(bk_comp_matrix), Intent(InOut) :: a
      Type(bk_comp_matrix), Intent(InOut) :: b
      ! Local Variables
      Integer :: ierr
      If(.not.Allocated(a%data_d)) Then
        Allocate(a%data_d(b%n,b%m),stat=ierr)
        cuda_check(ierr)
      Endif
      a%data_d = b%data_d
      Deallocate(b%data_d,stat=ierr)
      cuda_check(ierr)
      b%n = INVALID
      b%m = INVALID
      print_device_memory('bk_comp_matrix_copy')
    End Subroutine bk_comp_matrix_copy

    Subroutine bk_real_sort_eval_evecs(evals,evecs)
    Use sort
      Type(bk_vector_1D)  , Intent(InOut)           :: evals
      Type(bk_real_matrix), Intent(InOut), Optional :: evecs
      ! Local Variables
      Type(dim3) :: grid_size, block_size
      Integer    :: ldn
      Integer    :: ierr
      ! Device Variables
      Integer   , Dimension(:)  , Allocatable, MANAGED :: indices
      Real(fp64), Dimension(:,:), Allocatable, Device  :: mat_tmp
      ldn = Size(evals%data_d)
      Allocate(indices(ldn),stat=ierr)
      cuda_check(ierr)
      Call fsort(evals%data_d,indices,ldn,.True.)
      If(Present(evecs)) Then
        Allocate(mat_tmp(ldn,ldn),stat=ierr)
        cuda_check(ierr)
        mat_tmp = evecs%data_d
        block_size = dim3( MAX_THREADS_2D, MAX_THREADS_2D, 1 )
        grid_size = dim3( (ldn+MAX_THREADS_2D-1)/MAX_THREADS_2D, (ldn+MAX_THREADS_2D-1)/MAX_THREADS_2D, 1 )
        Call reorder_real_kernel<<<grid_size,block_size>>>(ldn,mat_tmp,evecs%data_d,indices)
        Deallocate(mat_tmp,stat=ierr)
        cuda_check(ierr)
      Endif
    End Subroutine bk_real_sort_eval_evecs

    Subroutine bk_comp_sort_eval_evecs(evals,evecs)
    Use sort
      Type(bk_vector_1D)  , Intent(InOut) :: evals
      Type(bk_comp_matrix), Intent(InOut) :: evecs
      ! Local Variables
      Type(dim3) :: grid_size, block_size
      Integer    :: ldn
      Integer    :: ierr
      ! Device Variables
      Integer      , Dimension(:)  , Allocatable, MANAGED :: indices
      Complex(fp64), Dimension(:,:), Allocatable, Device  :: mat_tmp
      ldn = Size(evals%data_d)
      Allocate(indices(ldn),stat=ierr)
      cuda_check(ierr)
      Call fsort(evals%data_d,indices,ldn,.True.)
      Allocate(mat_tmp(ldn,ldn),stat=ierr)
      cuda_check(ierr)
      mat_tmp = evecs%data_d
      block_size = dim3( MAX_THREADS_2D, MAX_THREADS_2D, 1 )
      grid_size = dim3( (ldn+MAX_THREADS_2D-1)/MAX_THREADS_2D, (ldn+MAX_THREADS_2D-1)/MAX_THREADS_2D, 1 )
      Call reorder_comp_kernel<<<grid_size,block_size>>>(ldn,mat_tmp,evecs%data_d,indices)
      Deallocate(mat_tmp,stat=ierr)
      cuda_check(ierr)
    End Subroutine bk_comp_sort_eval_evecs




    !! =========================================================================
    !! ft_data type bounded procedures - START
    !! =========================================================================
    Subroutine ft_data_set( this_map,k_grid,k_grid_mult,k_grid_lcm,n_couples, &
                             & sh_ptr,sh1_ptr,sh2_ptr,sh_sizes,pg_sh_ptr,        &
                             & stars_per_c,g_per_star,g_ptr,g_labels,g_vectors ) 
      Class(ft_data),                 Intent(InOut) :: this_map
      Integer       , Dimension(3)  , Intent(In   ) :: k_grid
      Integer       , Dimension(3)  , Intent(In   ) :: k_grid_mult
      Integer       ,                 Intent(In   ) :: k_grid_lcm
      Integer       ,                 Intent(In   ) :: n_couples
      Integer       , Dimension(:)  , Intent(In   ) :: sh_ptr
      Integer       , Dimension(:)  , Intent(In   ) :: sh1_ptr
      Integer       , Dimension(:)  , Intent(In   ) :: sh2_ptr
      Integer       , Dimension(:)  , Intent(In   ) :: sh_sizes
      Integer       , Dimension(:)  , Intent(In   ) :: pg_sh_ptr
      Integer       , Dimension(:)  , Intent(In   ) :: stars_per_c
      Integer       , Dimension(:)  , Intent(In   ) :: g_per_star
      Integer       , Dimension(:)  , Intent(In   ) :: g_ptr
      Integer       , Dimension(:)  , Intent(In   ) :: g_labels
      Integer       , Dimension(:,:), Intent(In   ) :: g_vectors
      ! Local Variables
      Integer :: couple, op
      Integer :: sh1, sh2, sh1_symm, sh2_symm
      Integer :: star, g
      Integer :: nugir
      Integer :: recursions, max_recursions = 0
      Integer :: start, end
      Integer :: ierr

      If(this_map%status.gt.STATUS_INIT) Then
        Call this_map%free()
      Endif

      this_map%k_grid      = k_grid       !! Shrinking factors
      this_map%k_grid_mult = k_grid_mult  !! Shrinking factors multipliers (to get MCM)
      this_map%k_grid_lcm  = k_grid_lcm   !! Least common multiplier of shrinking factors (MCM)
      this_map%n_couples   = n_couples    !! number of selected couples
      Allocate( this_map%shell1_size   (n_couples),stat=ierr )
      Allocate( this_map%shell1_start  (n_couples),stat=ierr )
      Allocate( this_map%shell2_size   (n_couples),stat=ierr )
      Allocate( this_map%shell2_start  (n_couples),stat=ierr )
      Allocate( this_map%pg_shell_start(n_couples),stat=ierr )
      Allocate( this_map%repetitions   (n_couples),stat=ierr )
      Do couple = 1, n_couples
        sh1 = sh1_ptr( sh_ptr(couple) + 1 )                  !! <---------- SHELL 1
        this_map%shell1_size (couple) = sh_sizes(sh1)
        this_map%shell1_start(couple) = Sum( sh_sizes( 1 : sh1-1 ))
        sh2 = sh2_ptr( sh_ptr(couple) + 1 )                  !! <---------- SHELL 2
        this_map%shell2_size (couple) = sh_sizes(sh2)
        this_map%shell2_start(couple) = Sum( sh_sizes( 1 : sh2-1 ))
        this_map%pg_shell_start(couple) = pg_sh_ptr(couple)  !! <---------- PGirr START SHELL POSITION
        recursions = 0                                       !! <---------- g VECS FOR EACH CONTRIBUTION
        Do star = 1, stars_per_c(couple)
          Do g = 1, g_per_star( g_ptr(couple) + star )
            recursions = recursions + 1
          Enddo
        Enddo
        this_map%repetitions(couple) = recursions
      Enddo ! shell loop
      max_recursions = Maxval(this_map%repetitions)
      Allocate( this_map%g_vecs(3,max_recursions,n_couples) )
      this_map%g_vecs = -Huge(INVALID)
      Do couple = 1, n_couples
        nugir = g_ptr(couple)
        recursions = 1
        Do star = 1, stars_per_c(couple)
          Do g = 1, g_per_star( g_ptr(couple) + star )
            nugir = nugir + 1
            this_map%g_vecs(1,recursions,couple) = g_vectors( 1, g_labels(nugir) )
            this_map%g_vecs(2,recursions,couple) = g_vectors( 2, g_labels(nugir) )
            this_map%g_vecs(3,recursions,couple) = g_vectors( 3, g_labels(nugir) )
            recursions = recursions + 1
          Enddo
        Enddo
      Enddo ! shell loop
      this_map%status = STATUS_DEVICE
    End Subroutine ft_data_set

    Subroutine ft_data_free(this_map)
      Class(ft_data), Intent(InOut) :: this_map
      this_map%k_grid      = [INVALID,INVALID,INVALID]
      this_map%k_grid_mult = [INVALID,INVALID,INVALID]
      this_map%k_grid_lcm  = INVALID
      this_map%n_couples   = INVALID
      If( Allocated(this_map%shell1_size   ) ) Deallocate(this_map%shell1_size   )
      If( Allocated(this_map%shell1_start  ) ) Deallocate(this_map%shell1_start  )
      If( Allocated(this_map%shell2_size   ) ) Deallocate(this_map%shell2_size   )
      If( Allocated(this_map%shell2_start  ) ) Deallocate(this_map%shell2_start  )
      If( Allocated(this_map%pg_shell_start) ) Deallocate(this_map%pg_shell_start)
      If( Allocated(this_map%repetitions   ) ) Deallocate(this_map%repetitions   )
      If( Allocated(this_map%g_vecs        ) ) Deallocate(this_map%g_vecs        )
      this_map%status = STATUS_INIT
    End Subroutine ft_data_free

    Subroutine ft_data_print(this_map,cry_memory,out)
      Class(ft_data) , Intent(In) :: this_map
      Integer        , Intent(In) :: cry_memory
      Integer        , Intent(In) :: out
      Write(out,"(1X,T26,A)") "*** FT DATA MAP MEMORY USAGE ***"
      Write(out,"(1X,'N SHELL COUPLES',T21,':   ',I10,T41,'MAX SHELL SIZE:',T69,I10)") &
           &    this_map%n_couples, Max( Maxval(this_map%shell1_size), Maxval(this_map%shell2_size) )
      Write(out,"(1X,'MAX REPETITIONS',T21,':   ',I10)") &
           &    Size(this_map%g_vecs, Dim=2)
      Write(out,"(1X,'CRYSTAL POINTERS   :   ',F10.2,' MB   (just as reference)')") Real(cry_memory,kind=fp64) * 4._fp64 / 1024._fp64 / 1024._fp64
      Write(out,"(1X,'NEW DEVICE DATA MAP:   ',F10.2,' MB')") Real( this_map%size(), kind=fp64 ) / 1024._fp64 / 1024._fp64 
      Write(out,"(1X,79('*'))")
    End Subroutine ft_data_print

    Function ft_data_get_alloc_mem(this_map) Result(mem)
    !! Returns the current allocated memory (integer) of self in bytes
      Class(ft_data), Intent(In   ) :: this_map
      Integer(int64), Intent(  Out) :: mem
      mem =                                              Size(this_map%k_grid        ) * Storage_size(this_map%k_grid        )
      mem =                                        mem + Size(this_map%k_grid_mult   ) * Storage_size(this_map%k_grid_mult   )
      mem =                                        mem +                                 Storage_size(this_map%k_grid_lcm    )
      mem =                                        mem +                                 Storage_size(this_map%n_couples     )
      If(Allocated(this_map%shell1_size   )) mem = mem + Size(this_map%shell1_size   ) * Storage_size(this_map%shell1_size   )
      If(Allocated(this_map%shell1_start  )) mem = mem + Size(this_map%shell1_start  ) * Storage_size(this_map%shell1_start  )
      If(Allocated(this_map%shell2_size   )) mem = mem + Size(this_map%shell2_size   ) * Storage_size(this_map%shell2_size   )
      If(Allocated(this_map%shell2_start  )) mem = mem + Size(this_map%shell2_start  ) * Storage_size(this_map%shell2_start  )
      If(Allocated(this_map%pg_shell_start)) mem = mem + Size(this_map%pg_shell_start) * Storage_size(this_map%pg_shell_start)
      If(Allocated(this_map%repetitions   )) mem = mem + Size(this_map%repetitions   ) * Storage_size(this_map%repetitions   )
      If(Allocated(this_map%g_vecs        )) mem = mem + Size(this_map%g_vecs        ) * Storage_size(this_map%g_vecs        )
      mem = mem / 8_int64
    End Function ft_data_get_alloc_mem
    !! =========================================================================
    !! ft_data type bounded procedures - END
    !! =========================================================================

    Subroutine bk_real_matrix_compute_density( c,nocc,pg,weights,coords )
      Type(bk_real_matrix),               Intent(InOut) :: c
      Integer             ,               Intent(In   ) :: nocc
      Real(fp64)          , Dimension(:), Intent(InOut) :: pg
      Real(fp64)          , Dimension(:), Intent(In   ) :: weights
      Integer             , Dimension(3), Intent(In   ) :: coords
      ! Local Variables
      Integer, Dimension(3) :: k_vec
      Integer               :: ierr
      Type(dim3)            :: block_size
      ! Device Variables
      Real(fp64), Dimension(:), Allocatable, Device :: weights_d
      Real(fp64), Dimension(:), Allocatable, Device :: density_d
      Integer   , Dimension(3),              Device :: kvec_d
      k_vec(1:3) = coords(1:3) * ft_map%k_grid_mult(1:3)
      Allocate(weights_d(nocc),stat=ierr)
      cuda_check(ierr)
      Allocate(density_d(Size(pg)))
      cuda_check(ierr)
      weights_d = weights
      density_d = pg
      kvec_d = k_vec
      block_size = dim3(Maxval(ft_map%shell1_size),Maxval(ft_map%shell2_size),1)
      If(block_size%x*block_size%y.gt.1024) print*,'ERRRRRRRRRRRRRRRROOOOOOOOOORRRRRRRRRR'
      call __CUDA_nvtx_profiler_start("compute_density_real",1)
      Call real_density_contribution_kernel<<<ft_map%n_couples,block_size>>> &
           & ( c%data_d,nocc,ft_map%shell1_size,ft_map%shell2_size, &
           &   ft_map%shell1_start,ft_map%shell2_start,weights_d,density_d, &
           &   ft_map%pg_shell_start,ft_map%repetitions,kvec_d,ft_map%k_grid_lcm,ft_map%g_vecs )
      cuda_get_error_check
      call __CUDA_nvtx_profiler_end
      pg = density_d  
      Deallocate(weights_d,density_d,stat=ierr)
      cuda_check(ierr)
    End Subroutine bk_real_matrix_compute_density

    Subroutine bk_comp_matrix_compute_density( c,nocc,pg,weights,coords )
      Type(bk_comp_matrix),               Intent(InOut) :: c
      Integer             ,               Intent(In   ) :: nocc
      Real(fp64)          , Dimension(:), Intent(InOut) :: pg
      Real(fp64)          , Dimension(:), Intent(In   ) :: weights
      Integer             , Dimension(3), Intent(In   ) :: coords
      ! Local Variables
      Integer, Dimension(3) :: k_vec
      Integer               :: ierr
      Type(dim3)            :: block_size
      ! Device Variables
      Real(fp64), Dimension(:), Allocatable, Device :: weights_d
      Real(fp64), Dimension(:), Allocatable, Device :: density_d
      Integer   , Dimension(3),              Device :: kvec_d
      k_vec(1:3) = coords(1:3) * ft_map%k_grid_mult(1:3)
      Allocate(weights_d(nocc),stat=ierr)
      cuda_check(ierr)
      Allocate(density_d(Size(pg)))
      cuda_check(ierr)
      weights_d = weights
      density_d = pg
      kvec_d = k_vec
      block_size = dim3(Maxval(ft_map%shell1_size),Maxval(ft_map%shell2_size),1)
      If(block_size%x*block_size%y.gt.1024) print*,'ERRRRRRRRRRRRRRRROOOOOOOOOORRRRRRRRRR'
      call __CUDA_nvtx_profiler_start("compute_density_comp",1)
      Call comp_density_contribution_kernel<<<ft_map%n_couples,block_size>>> &
           & ( c%data_d,nocc,ft_map%shell1_size,ft_map%shell2_size, &
           &   ft_map%shell1_start,ft_map%shell2_start,weights_d,density_d, &
           &   ft_map%pg_shell_start,ft_map%repetitions,kvec_d,ft_map%k_grid_lcm,ft_map%g_vecs )
      cuda_get_error_check
      call __CUDA_nvtx_profiler_end
      pg = density_d  
      Deallocate(weights_d,density_d,stat=ierr)
      cuda_check(ierr)
    End Subroutine bk_comp_matrix_compute_density












!! --- Matrix Related ---
    Function bk_real_allocated(matrix) Result(res)
      Logical                          :: res
      Type(bk_real_matrix), Intent(In) :: matrix
      res = Allocated(matrix%data_d)
    End Function bk_real_allocated

    Function bk_comp_allocated(matrix) Result(res)
      Logical                          :: res
      Type(bk_comp_matrix), Intent(In) :: matrix
      res = Allocated(matrix%data_d)
    End Function bk_comp_allocated

    Subroutine bk_real_shift_diag(matrix,shift,nocc)
      Type(bk_real_matrix), Intent(InOut) :: matrix
      Real(fp64)          , Intent(In   ) :: shift
      Integer             , Intent(In   ) :: nocc
      ! Local Variables
      Integer :: ldm
      Integer :: m
      Integer :: i, j
      ! Device Variables
      Real(fp64), Dimension(:,:), Allocatable, Device :: temp
      ldm = matrix%n
      Allocate(temp(ldm,ldm))
      temp = matrix%data_d
      !$cuf kernel do(1) <<< *, * >>>
      Do i = 1, nocc
        temp(i,i) = temp(i,i) - shift
      Enddo
      matrix%data_d = temp
    End Subroutine bk_real_shift_diag

    Subroutine bk_comp_shift_diag(matrix,shift,nocc)
    !! Just shift the real component on the diagonal
      Type(bk_comp_matrix), Intent(InOut) :: matrix
      Real(fp64)          , Intent(In   ) :: shift
      Integer             , Intent(In   ) :: nocc
      ! Local Variables
      Integer :: ldm
      Integer :: m
      Integer :: i, j
      ! Device Variables
      Complex(fp64), Dimension(:,:), Allocatable, Device :: temp
      ldm = matrix%n
      Allocate(temp(ldm,ldm))
      temp = matrix%data_d
      !$cuf kernel do(1) <<< *, * >>>
      Do i = 1, nocc
        temp(i,i) = temp(i,i) - Cmplx(shift,0.0_fp64,Kind=fp64)
      Enddo
      matrix%data_d = temp
    End Subroutine bk_comp_shift_diag

    Subroutine bk_compute_real_tester_sub(matrix,nocc,tstr)
      Type(bk_real_matrix), Intent(In   ) :: matrix
      Integer             , Intent(In   ) :: nocc
      Real(fp64)          , Intent(  Out) :: tstr
      ! Local Variables
      Integer    :: ldm
      Type(dim3) :: grid_size, block_size
      ! Device Variables
      Real(fp64), Device :: tstr_d
      ldm = Max(matrix%n,matrix%m)
      tstr_d = 0.0_fp64
      block_size = dim3(MAX_THREADS_2D,MAX_THREADS_2D,1)
      grid_size = dim3( (ldm-nocc+MAX_THREADS_2D-1)/MAX_THREADS_2D, & ! x: nocc+1:ldm
                      & (nocc+MAX_THREADS_2D-1)/MAX_THREADS_2D, 1 )   ! y: 1:nocc
      Call real_compute_tester<<<grid_size,block_size>>>(matrix%data_d,nocc,ldm,tstr_d)
      tstr = tstr_d
    End Subroutine bk_compute_real_tester_sub

    Subroutine bk_compute_complex_tester_sub(matrix,nocc,tstr)
    !! Just computing tester value on the real part of complex matrices
      Type(bk_comp_matrix), Intent(In   ) :: matrix
      Integer             , Intent(In   ) :: nocc
      Real(fp64)          , Intent(  Out) :: tstr
      ! Local Variables
      Integer    :: ldm
      Type(dim3) :: grid_size, block_size
      ! Device Variables
      Real(fp64), Device :: tstr_d
      ldm = Max(matrix%n,matrix%m)
      tstr_d = 0.0_fp64
      block_size = dim3(MAX_THREADS_2D,MAX_THREADS_2D,1)
      grid_size = dim3( (ldm-nocc+MAX_THREADS_2D-1)/MAX_THREADS_2D, & ! x: nocc+1:ldm
                      & (nocc+MAX_THREADS_2D-1)/MAX_THREADS_2D, 1 )   ! y: 1:nocc
      Call comp_compute_tester<<<grid_size,block_size>>>(matrix%data_d,nocc,ldm,tstr_d)
      tstr = tstr_d
    End Subroutine bk_compute_complex_tester_sub

    Subroutine bk_shift_1d(vec,nocc,shift)
      Type(bk_vector_1D), Intent(InOut) :: vec
      Integer           , Intent(In   ) :: nocc
      Real(fp64)        , Intent(In   ) :: shift
      ! Local Variables 
      Integer    :: ierr
      Type(dim3) :: grid_size,block_size
      block_size = dim3( MAX_THREADS, 1, 1 )
      grid_size  = dim3( ((nocc + MAX_THREADS - 1) / MAX_THREADS), 1, 1 )
      Call vector_shift_kernel<<<grid_size,block_size>>>(vec%data_d,nocc,shift)
    End Subroutine bk_shift_1d

    Subroutine bk_host_real_alloc(a,n,m)
      Type(bk_real_matrix), Intent(InOut) :: a
      Integer             , Intent(In   ) :: n
      Integer             , Intent(In   ) :: m
      ! Local Variables
      Integer :: ierr
      a%n = n
      a%m = m
      Allocate(a%data(n,m), stat=ierr)
    End Subroutine bk_host_real_alloc

    Subroutine bk_host_comp_alloc(a,n,m)
      Type(bk_comp_matrix), Intent(InOut) :: a
      Integer             , Intent(In   ) :: n
      Integer             , Intent(In   ) :: m
      ! Local Variables
      Integer :: ierr
      a%n = n
      a%m = m
      Allocate(a%data(n,m), stat=ierr)
    End Subroutine bk_host_comp_alloc

    Subroutine bk_host_real_dealloc(a)
      Type(bk_real_matrix), Intent(InOut) :: a
      ! Local Variables
      Integer :: ierr
      If(.not.Allocated(a%data)) cuda_error("deallocating not-allocated-variable")
      Deallocate(a%data, stat=ierr)
      a%n = INVALID
      a%m = INVALID
    End Subroutine bk_host_real_dealloc

    Subroutine bk_host_comp_dealloc(a)
      Type(bk_comp_matrix), Intent(InOut) :: a
      ! Local Variables
      Integer :: ierr
      If(.not.Allocated(a%data)) cuda_error("deallocating not-allocated-variable")
      Deallocate(a%data, stat=ierr)
      a%n = INVALID
      a%m = INVALID
    End Subroutine bk_host_comp_dealloc

    Subroutine bk_host_real_matrix_copy(a,b)
      Type(bk_real_matrix), Intent(InOut) :: a
      Type(bk_real_matrix), Intent(InOut) :: b
      ! Local Variables
      Integer :: ierr
      If(.not.Allocated(a%data)) Allocate(a%data(b%n,b%m),stat=ierr)
      a%data = b%data
      Deallocate(b%data,stat=ierr)
      b%n = INVALID
      b%m = INVALID
    End Subroutine bk_host_real_matrix_copy

    Subroutine bk_host_comp_matrix_copy(a,b)
      Type(bk_comp_matrix), Intent(InOut) :: a
      Type(bk_comp_matrix), Intent(InOut) :: b
      ! Local Variables
      Integer :: ierr
      If(.not.Allocated(a%data)) Allocate(a%data(b%n,b%m),stat=ierr)
      a%data = b%data
      Deallocate(b%data,stat=ierr)
      b%n = INVALID
      b%m = INVALID
    End Subroutine bk_host_comp_matrix_copy

    Subroutine bk_host_set_raw_real(a,raw_data)
      Type(bk_real_matrix),                 Intent(InOut) :: a
      Real(fp64)          , Dimension(:,:), Intent(In   ) :: raw_data
      If(Size(raw_data,Dim=1).ne.a%n) cuda_error("wrong domension")
      If(Size(raw_data,Dim=2).ne.a%m) cuda_error("wrong dimension")
      a%data = raw_data
    End Subroutine bk_host_set_raw_real

    Subroutine bk_host_set_raw_comp(a,raw_data)
      Type(bk_comp_matrix),                 Intent(InOut) :: a
      Complex(fp64)       , Dimension(:,:), Intent(In   ) :: raw_data
      If(Size(raw_data,Dim=1).ne.a%n) cuda_error("wrong domension")
      If(Size(raw_data,Dim=2).ne.a%m) cuda_error("wrong dimension")
      a%data = raw_data
    End Subroutine bk_host_set_raw_comp

    Subroutine bk_host_get_raw_real(a,raw_data)
      Type(bk_real_matrix),                              Intent(InOut) :: a
      Real(fp64)          , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      ! Local Variables
      Integer :: ldn, ldm
      ldn = a%n
      ldm = a%m
      Allocate(raw_data(ldn,ldm))
      raw_data = a%data
    End Subroutine bk_host_get_raw_real

    Subroutine bk_host_get_raw_comp(a,raw_data)
      Type(bk_comp_matrix),                              Intent(InOut) :: a
      Complex(fp64)       , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      ! Local Variables
      Integer :: ldn, ldm
      ldn = a%n
      ldm = a%m
      Allocate(raw_data(ldn,ldm))
      raw_data = a%data
    End Subroutine bk_host_get_raw_comp

    Subroutine bk_real_host_to_dev(mat)
      Type(bk_real_matrix), Intent(InOut) :: mat
      mat%data_d = mat%data
    End Subroutine bk_real_host_to_dev

    Subroutine bk_comp_host_to_dev(mat)
      Type(bk_comp_matrix), Intent(InOut) :: mat
      mat%data_d = mat%data
    End Subroutine bk_comp_host_to_dev

    Subroutine bk_real_dev_to_host(mat)
      Type(bk_real_matrix), Intent(InOut) :: mat
      mat%data = mat%data_d
    End Subroutine bk_real_dev_to_host

    Subroutine bk_comp_dev_to_host(mat)
      Type(bk_comp_matrix), Intent(InOut) :: mat
      mat%data = mat%data_d
    End Subroutine bk_comp_dev_to_host






































!! -----------------------------------------------------------------------
!!  Private Procedures
!! -----------------------------------------------------------------------

    Subroutine __CUDA_safe_call(ierr,name,line)            !! cuda_safe_Call
      Integer         , Intent(In) :: ierr
      Character(Len=*), Intent(In) :: name
      Integer         , Intent(In) :: line
      #ifdef DASH_BACK_DB
      If(ierr.ne.cudaSuccess) Then
        Write(*,"(A,I0,': ',A)") " ***** CUDA ERROR ***** in file "//trim(name)//" at line ",&
        & line, cudaGetErrorString(ierr)
      Endif
      #endif
    End Subroutine __CUDA_safe_call

    Subroutine __CUDA_error_msg(msg,name,line)             !! cuda_error
      Character(Len=*), Intent(In) :: msg
      Character(Len=*), Intent(In) :: name
      Integer         , Intent(In) :: line
      Write(*,"(A,I0,': ',A)") " ***** CUDA ERROR ***** in file "//trim(name)//" at line ", line, msg
    End Subroutine __CUDA_error_msg

    Subroutine __CUDA_error_check(ierr,name,line)          !! cuda_check
      Integer         , Intent(In) :: ierr
      Character(Len=*), Intent(In) :: name
      Integer         , Intent(In) :: line
      #ifdef DASH_BACK_DB
      If(ierr.ne.cudaSuccess) Then
        Write(*,"(A,I0,': ',A)") " ***** CUDA ERROR ***** in file "//trim(name)//" at line ",&
        & line, cudaGetErrorString(ierr)
      Endif
      #endif
    End Subroutine __CUDA_error_check

    Subroutine __CUDA_lib_error_check(ierr,name,line)      !! culib_check
      Integer         , Intent(In) :: ierr
      Character(Len=*), Intent(In) :: name
      Integer         , Intent(In) :: line
      #ifdef DASH_BACK_DB
      If(ierr.ne.cudaSuccess) Then
        Write(*,"(A,I0,': ',I0)") " ***** CUDA ERROR ***** in file "//trim(name)//" at line ",&
        & line, ierr
      Endif
      #endif
    End Subroutine __CUDA_lib_error_check

    Subroutine __CUDA_error_check_var(name,line)           !! cuda_get_error_check
      Character(Len=*), Intent(In) :: name
      Integer         , Intent(In) :: line
      ! Local Variables
      Integer :: ierr
      #ifdef DASH_BACK_DB
      ierr = cudaGetLastError()
      If(ierr.ne.cudaSuccess) Then
        Write(*,"(A,I0,': ',A)") " ***** CUDA ERROR ***** in file "//trim(name)//" at line ",&
        & line, cudaGetErrorString(ierr)
      Endif
      #endif
    End Subroutine __CUDA_error_check_var

    Logical Function query_check()
    !! Check if device is assigned, if not DO NOT QUERY GPU...
      query_check = .False.
      If(assigned_dev.ne.INVALID) query_check = .True.
    End Function query_check



    Subroutine __CUDA_nvtx_profiler_start(label,i)
    #ifdef DASH_BACK_PROF
    Use nvtx
    #endif
      Character(Len=*), Intent(In)           :: label
      Integer         , Intent(In), Optional :: i
      ! Local Variables
      Integer :: tag = 16
      #ifdef DASH_BACK_PROF
      If(Present(i)) tag = i
      Call nvtxStartRange(label,tag)
      #endif
    End Subroutine __CUDA_nvtx_profiler_start

    Subroutine __CUDA_nvtx_profiler_end()
    #ifdef DASH_BACK_PROF
    Use nvtx
      Call nvtxEndRange()
    #endif
    End Subroutine __CUDA_nvtx_profiler_end




    Attributes(Global) Subroutine real_set_lower_kernel_v2(data,val,ldn,ldm)
      Real(fp64), Dimension(:,:), Intent(InOut) :: data
      Real(fp64), Value         , Intent(In   ) :: val
      Integer   , Value         , Intent(In   ) :: ldn
      Integer   , Value         , Intent(In   ) :: ldm
      ! Local Variables
      Integer :: row, col
      If(blockIdx%y.le.blockIdx%x) Then
        row = (blockIdx%x-1) * blockDim%x + threadIdx%x
        col = (blockIdx%y-1) * blockDim%y + threadIdx%y
        If( row.le.ldn .and. col.le.ldm ) Then
          If( row .gt. col ) data(row,col) = val
        Endif
      Endif
    End Subroutine real_set_lower_kernel_v2

    Attributes(Global) Subroutine comp_set_lower_kernel(data,val,ldn,ldm)
      Complex(fp64), Dimension(:,:), Intent(InOut) :: data
      Complex(fp64), Value         , Intent(In   ) :: val
      Integer      , Value         , Intent(In   ) :: ldn
      Integer      , Value         , Intent(In   ) :: ldm
      ! Local Variables
      Integer :: row, col
      If(blockIdx%y.le.blockIdx%x) Then
        row = (blockIdx%x-1) * blockDim%x + threadIdx%x
        col = (blockIdx%y-1) * blockDim%y + threadIdx%y
        If( row.le.ldn .and. col.le.ldm ) Then
          If( row .gt. col ) data(row,col) = val
        Endif
      Endif
    End Subroutine comp_set_lower_kernel

    Attributes(Global) Subroutine real_compute_tester(data,nocc,ldm,tstr)
      Real(fp64), Dimension(:,:), Intent(In   ) :: data
      Integer   , Value         , Intent(In   ) :: nocc
      Integer   , Value         , Intent(In   ) :: ldm
      Real(fp64), Device        , Intent(InOut) :: tstr
      ! Local Variables
      Integer    :: row, col
      Real(fp64) :: val, waste
      row = nocc + (blockIdx%x-1) * blockDim%x + threadIdx%x
      col = (blockIdx%y-1) * blockDim%y + threadIdx%y
      If( row.le.ldm .and. col.le.nocc ) Then
        val = data(row,col)**2 / ( Abs(data(row,row)-data(col,col))+0.004_fp64 )
        waste = atomicMax(tstr, val)
      Endif
    End Subroutine real_compute_tester

    Attributes(Global) Subroutine comp_compute_tester(data,nocc,ldm,tstr)
      Complex(fp64), Dimension(:,:), Intent(In   ) :: data
      Integer      , Value         , Intent(In   ) :: nocc
      Integer      , Value         , Intent(In   ) :: ldm
      Real(fp64)   , Device        , Intent(InOut) :: tstr
      ! Local Variables
      Integer      :: row, col
      Real(fp64) :: val, waste
      row = nocc + (blockIdx%x-1) * blockDim%x + threadIdx%x
      col = (blockIdx%y-1) * blockDim%y + threadIdx%y
      If( row.le.ldm .and. col.le.nocc ) Then
        val = Real(data(row,col),kind=fp64)**2 / &
            & ( Abs(Real(data(row,row),kind=fp64)-Real(data(col,col),kind=fp64))+0.004_fp64 )
        waste = atomicMax(tstr, val)
      Endif
    End Subroutine comp_compute_tester

    Attributes(Global) Subroutine vector_shift_kernel(data,nocc,shift)
      Real(fp64), Dimension(:), Intent(InOut) :: data
      Integer   , Value       , Intent(In   ) :: nocc
      Real(fp64), Value       , Intent(In   ) :: shift
      ! Local Variables
      Integer :: idx
      idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
      If(idx.le.nocc) data(idx) = data(idx) + shift
    End Subroutine vector_shift_kernel

    Attributes(Global) Subroutine get_bands_range_kernel(levels,bmax_val,bmin_val,bmax_k,bmin_k)
      Real(fp64), Dimension(:,:,:), Intent(In   ) :: levels
      Real(fp64), Dimension(:,:)  , Intent(InOut) :: bmax_val
      Real(fp64), Dimension(:,:)  , Intent(InOut) :: bmin_val
      Integer   , Dimension(:,:)  , Intent(InOut) :: bmax_k
      Integer   , Dimension(:,:)  , Intent(InOut) :: bmin_k
      ! Local Variables
      Integer :: band_idx
      Integer :: spin_idx
      Integer :: k_idx

      band_idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
      spin_idx = blockIdx%y

      If( band_idx.le.Size(levels,Dim=1) ) Then
        Do k_idx = 1, Size(levels,Dim=3)
          If( levels(band_idx,spin_idx,k_idx) .gt. bmax_val(band_idx,spin_idx) ) Then !! Upper limit
            bmax_val(band_idx,spin_idx) = levels(band_idx,spin_idx,k_idx)
            bmax_k(band_idx,spin_idx)   = k_idx
          Endif
          If( levels(band_idx,spin_idx,k_idx) .lt. bmin_val(band_idx,spin_idx) ) Then !! Lower limit
            bmin_val(band_idx,spin_idx) = levels(band_idx,spin_idx,k_idx)
            bmin_k(band_idx,spin_idx)   = k_idx
          Endif
        Enddo
      Endif
    End Subroutine get_bands_range_kernel

    Attributes(Global) Subroutine reorder_real_kernel(ldn,mat_in,mat_out,order)
      Integer   , Value         , Intent(In   ) :: ldn
      Real(fp64), Dimension(:,:), Intent(In   ) :: mat_in
      Real(fp64), Dimension(:,:), Intent(InOut) :: mat_out
      Integer   , Dimension(:)  , Intent(In   ) :: order
      ! Local Variables
      Integer :: row, col, old_col
      row = threadIdx%x + (blockIdx%x - 1) * blockDim%x
      col = threadIdx%y + (blockIdx%y - 1) * blockDim%y
      If(row.le.ldn .and. col.le.ldn) Then
        old_col = order(col)
        mat_out(row,col) = mat_in(row,old_col)
      Endif
    End Subroutine reorder_real_kernel

    Attributes(Global) Subroutine reorder_comp_kernel(ldn,mat_in,mat_out,order)
      Integer      , Value         , Intent(In   ) :: ldn
      Complex(fp64), Dimension(:,:), Intent(In   ) :: mat_in
      Complex(fp64), Dimension(:,:), Intent(InOut) :: mat_out
      Integer      , Dimension(:)  , Intent(In   ) :: order
      ! Local Variables
      Integer :: row, col, old_col
      row = threadIdx%x + (blockIdx%x - 1) * blockDim%x
      col = threadIdx%y + (blockIdx%y - 1) * blockDim%y
      If(row.le.ldn .and. col.le.ldn) Then
        old_col = order(col)
        mat_out(row,col) = mat_in(row,old_col)
      Endif
    End Subroutine reorder_comp_kernel


    Attributes(Global) Subroutine real_density_contribution_kernel &
      & ( data,nocc,sh1,sh2,sh1_st,sh2_st,weights,dens,pg_st,reps,k_vec,lcm,g_vec )
      Real(fp64), Dimension(:,:)  , Intent(InOut) :: data
      Integer   , Value           , Intent(In   ) :: nocc                 
      Integer   , Dimension(:)    , Intent(In   ) :: sh1
      Integer   , Dimension(:)    , Intent(In   ) :: sh2
      Integer   , Dimension(:)    , Intent(In   ) :: sh1_st
      Integer   , Dimension(:)    , Intent(In   ) :: sh2_st
      Real(fp64), Dimension(:)    , Intent(In   ) :: weights
      Real(fp64), Dimension(:)    , Intent(InOut) :: dens
      Integer   , Dimension(:)    , Intent(In   ) :: pg_st
      Integer   , Dimension(:)    , Intent(In   ) :: reps
      Integer   , Dimension(3)    , Intent(In   ) :: k_vec
      Integer   , Value           , Intent(In   ) :: lcm
      Integer   , Dimension(:,:,:), Intent(In   ) :: g_vec
      ! Local Variables
      Integer    :: size1
      Integer    :: size2
      Integer    :: start1
      Integer    :: start2
      Integer    :: pos
      Integer    :: r
      Real(fp64) :: expo
      Integer    :: k_dot_g
      Real(fp64) :: val

      size1 = sh1(blockIdx%x)
      size2 = sh2(blockIdx%x)
      start1 = sh1_st(blockIdx%x)
      start2 = sh2_st(blockIdx%x)

      If( threadIdx%x.le.size1 .and. threadIdx%y.le.size2 ) Then
        val = Sum( data(start1+threadIdx%x,1:nocc) * data(start2+threadIdx%y,1:nocc) * weights(1:nocc) )
        pos = pg_st(blockIdx%x) + threadIdx%y + ( (threadIdx%x-1)*size2 )
        Do r = 1, reps(blockIdx%x) 
          k_dot_g = k_vec(1)*g_vec(1,r,blockIdx%x) + k_vec(2)*g_vec(2,r,blockIdx%x) + k_vec(3)*g_vec(3,r,blockIdx%x)
          expo = Cos( Real(Modulo(k_dot_g,lcm),kind=fp64) * 2.0_fp64*pi/Real(lcm,kind=fp64) )
          dens(pos) = dens(pos) + ( val * expo )
          pos = pos + size1*size2
        Enddo
      Endif

    End Subroutine real_density_contribution_kernel

    Attributes(Global) Subroutine comp_density_contribution_kernel &
      & ( data,nocc,sh1,sh2,sh1_st,sh2_st,weights,dens,pg_st,reps,k_vec,lcm,g_vec )
      Complex(fp64), Dimension(:,:)  , Intent(InOut) :: data
      Integer      , Value           , Intent(In   ) :: nocc                 
      Integer      , Dimension(:)    , Intent(In   ) :: sh1
      Integer      , Dimension(:)    , Intent(In   ) :: sh2
      Integer      , Dimension(:)    , Intent(In   ) :: sh1_st
      Integer      , Dimension(:)    , Intent(In   ) :: sh2_st
      Real(fp64)   , Dimension(:)    , Intent(In   ) :: weights
      Real(fp64)   , Dimension(:)    , Intent(InOut) :: dens
      Integer      , Dimension(:)    , Intent(In   ) :: pg_st
      Integer      , Dimension(:)    , Intent(In   ) :: reps
      Integer      , Dimension(3)    , Intent(In   ) :: k_vec
      Integer      , Value           , Intent(In   ) :: lcm
      Integer      , Dimension(:,:,:), Intent(In   ) :: g_vec
      ! Local Variables
      Integer    :: size1
      Integer    :: size2
      Integer    :: start1
      Integer    :: start2
      Integer    :: pos
      Integer    :: r
      Real(fp64) :: expo_cos, expo_sen
      Integer    :: k_dot_g
      Real(fp64) :: val_r, val_i

      size1 = sh1(blockIdx%x)
      size2 = sh2(blockIdx%x)
      start1 = sh1_st(blockIdx%x)
      start2 = sh2_st(blockIdx%x)

      If( threadIdx%x.le.size1 .and. threadIdx%y.le.size2 ) Then

        !! 1r*2r + 1i*2i
        val_r = Sum(                                                &
              & ( Real(data(start1+threadIdx%x,1:nocc),kind=fp64) * &
              &   Real(data(start2+threadIdx%y,1:nocc),kind=fp64) * &
              &   weights(1:nocc) )                               + &
              & ( Aimag(data(start1+threadIdx%x,1:nocc))          * &
              &   Aimag(data(start2+threadIdx%y,1:nocc))          * &
              &   weights(1:nocc) )                                 )

        !! 1i*2r - 1r*2i
        val_i = Sum( &
              & ( Aimag(data(start1+threadIdx%x,1:nocc))          * &
              &   Real(data(start2+threadIdx%y,1:nocc),kind=fp64) * &
              &   weights(1:nocc) )                               - &
              & ( Real(data(start1+threadIdx%x,1:nocc),kind=fp64) * &
              &   Aimag(data(start2+threadIdx%y,1:nocc))          * &
              &   weights(1:nocc) )                                 )

        pos = pg_st(blockIdx%x) + threadIdx%y + ( (threadIdx%x-1)*size2 )

        Do r = 1, reps(blockIdx%x) 
          k_dot_g = k_vec(1)*g_vec(1,r,blockIdx%x) + k_vec(2)*g_vec(2,r,blockIdx%x) + k_vec(3)*g_vec(3,r,blockIdx%x)
          expo_cos = Cos( Real(Modulo(k_dot_g,lcm),kind=fp64) * 2.0_fp64*pi/Real(lcm,kind=fp64) )
          expo_sen = Sin( Real(Modulo(k_dot_g,lcm),kind=fp64) * 2.0_fp64*pi/Real(lcm,kind=fp64) )
          dens(pos) =  dens(pos) + (val_r * expo_cos) + (val_i * expo_sen)
          pos = pos + size1*size2
        Enddo
      Endif

    End Subroutine comp_density_contribution_kernel

  !! BASE CPU CODE FOR PDIG REAL
  !! ===========================================================================
    !! Do couple = 1, n_couples    !! LOOP ON SHELL COUPLES
    !!   sh1 = shells_1( shell_ptr(couple) + 1 )
    !!   sh2 = shells_2( shell_ptr(couple) + 1 )
    !!   ld_sh1 = shell_sizes(sh1)  !! size shell 1
    !!   ld_sh2 = shell_sizes(sh2)  !! size shell 2
    !!   ptr_sh1 = Sum( shell_sizes( 1 : sh1-1 ))
    !!   ptr_sh2 = Sum( shell_sizes( 1 : sh2-1 ))
    !!   Allocate( ev1( ld_sh1,occupied_states ) )
    !!   Allocate( ev2( ld_sh2,occupied_states ) )
    !!   Allocate( mat( ld_sh2,ld_sh1 ) ) !! already trsposed
    !!   ev1 = evecs( ptr_sh1+1:ptr_sh1+ld_sh1, 1:occupied_states )
    !!   ev2 = evecs( ptr_sh2+1:ptr_sh2+ld_sh2, 1:occupied_states )
    !!   do mu = 1, ld_sh2
    !!     ev2(mu,1:occupied_states) = ev2(mu,1:occupied_states) * weights(1:occupied_states)
    !!   enddo
    !!   !call dgemm('N', 'T', smu, snu, si, 1.0_fp64, ci, smu, cidag, snu, 0.0_fp64, my_a, smu)
    !!   !! To generate the transpose of my_a... should allocate my_a(snu,smu)
    !!   call dgemm('N', 'T', ld_sh2, ld_sh1, occupied_states, 1.0_fp64, ev2, ld_sh2, ev1, ld_sh1, 0.0_fp64, mat, ld_sh2)
    !!   Deallocate(ev1, ev2)
    !!   numbo = couple_start(couple) !+ nstap
    !!   nugir = g_vec_start(couple)
    !!   Do star = 1, g_stars_per_couple(couple)
    !!     Do g = 1, g_per_star( g_vec_start(couple) + star )
    !!          nugir = nugir + 1
    !!          mcm = this%k_grid_lcm
    !!          mg = NNGI(NUGIR)
    !!  !print*, cos( MOD(this_ks%coord(1) * LG(1,MG) + &
    !!  !              &  this_ks%coord(2) * LG(2,MG) + &
    !!  !              &  this_ks%coord(3) * LG(3,MG) + &
    !!  !              &  mcm*8192                       , mcm ) *  &
    !!  !              &  2.0_fp64*pi / mcm )
    !!       my_pg_irr( numbo+1 : numbo+ld_sh1*ld_sh2 ) = my_pg_irr( numbo+1 : numbo+ld_sh1*ld_sh2 ) &
    !!                               & + Reshape( mat(1:ld_sh2,1:ld_sh1), [ld_sh1*ld_sh2] )
    !!       !print*,numbo+1,'placing',mat
    !!       numbo = numbo + ld_sh1*ld_sh2
    !!     Enddo
    !!   Enddo
    !!   Deallocate(mat)
    !!   !if(couple.eq.4) exit
    !! Enddo
    !! Deallocate(evecs)
  !! ===========================================================================


















































    !Subroutine error_check(code,name,step)
    !!! - DEVICE QUERY -
    !!! Check if code is an error for cuda functions
    !  Integer         , Intent(In)           :: code
    !  Character(Len=*), Intent(In)           :: name
    !  Integer         , Intent(In), Optional :: step
    !  ! Local Variables
    !  If(code.ne.cudaSuccess) Then
    !    If(Present(step)) Then
    !      Write(*,"(A,I0,3X,'(',A,')')") &
    !      & "**CUDA** ERROR in "//trim(name)//" step: ",step,cudaGetErrorString(code)
    !    Else
    !      Write(*,"(A,3X,'(',A,')')") &
    !      & "**CUDA** ERROR in "//trim(name),cudaGetErrorString(code)
    !    Endif
    !  Endif
    !End Subroutine error_check

End Module backend_module
