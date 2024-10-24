Module backend_module

!===========================================================
!    @author: Giacomo Ambrogio       date: Sep 2024
!
!    DASH - Device Accelerated Solutions for HPC
!
!    Backend Module
!
!    This module is the only one which interact with GPUs.
!
!    This version is dedicated to CUDA Fortran and CUDA
!    libraries (NVIDIA GPUs only)
!
!===========================================================

Use base_numbers
Use cudafor
Use cublas_v2
Use cusolverdn

  Implicit None
  Private

  Logical, Private :: do_HOST_lapack = .False.
  Logical, Private :: do_DEV_cuda    = .True.

  !! --- VARIABLES ---
  Type(cublashandle)    , Private :: blas_handle
  Type(cusolverdnhandle), Private :: solver_handle
  Integer               , Private :: init_blas_handle_status   = 0
  Integer               , Private :: init_solver_handle_status = 0
  Type(cusolverdnparams), Private :: solver_params
  Integer               , Private :: init_solver_params_status   = 0

  !! --- DERIVED TYPES ---
  Type, Public :: bk_real_matrix
    Real(double), Dimension(:,:), Allocatable         :: data
    Real(double), Dimension(:,:), Device, Allocatable :: ddata
  End Type bk_real_matrix
  Type, Public :: bk_comp_matrix
    Complex(double), Dimension(:,:), Allocatable :: data
    Complex(double), Dimension(:,:), Device, Allocatable :: ddata
  End Type bk_comp_matrix
  Type, Public :: bk_vector_1D
    Real(double), Dimension(:), Allocatable :: data
    Real(double), Dimension(:), Allocatable, Device :: ddata
  End Type bk_vector_1D

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
  Public :: bk_set_raw_real     , bk_set_raw_comp
  Public :: bk_get_raw_real     , bk_get_raw_comp
  Public :: bk_get_raw_real_vect, bk_get_raw_comp_vect
  Public :: bk_set_raw_1d
  Public :: bk_get_raw_1d
  !! --- LA Operations ---
  Public :: bk_real_cholesky    , bk_comp_cholesky
  Public :: bk_real_set_lower_to, bk_comp_set_lower_to
  Public :: bk_real_invert      , bk_comp_invert
  Public :: bk_real_mult        , bk_comp_mult
  Public :: bk_real_sym_diag    , bk_comp_sym_diag
  Public :: bk_real_matrix_copy , bk_comp_matrix_copy
  !! --- Matrix Related ---
  Public :: bk_real_allocated         , bk_comp_allocated
  Public :: bk_real_shift_diag        , bk_comp_shift_diag
  Public :: bk_compute_real_tester_sub, bk_compute_complex_tester_sub

  Contains

!! --- Utility ---
    Integer Function bk_get_backend()
    !! return backend identifier index
      bk_get_backend = 2
    End Function bk_get_backend

    Integer Function bk_get_visible_devices(limit)
    !! return number of available decices for this process
      Integer, Intent(In), Optional :: limit
      ! Local Variables
      Integer :: ierr
      Integer :: n_dev
      ierr = cudaGetDeviceCount(n_dev)
      If(Present(limit))Then
        If(n_dev.gt.limit) n_dev = limit
      Endif
      bk_get_visible_devices = n_dev
    End Function bk_get_visible_devices

    Integer Function bk_get_current_device()
    !! return current device associated with this process
      ! Local Variables
      Integer :: ierr
      ierr =  cudaGetDevice(bk_get_current_device)
    End Function bk_get_current_device

    Subroutine bk_set_current_device(index)
    !! set gpu index for thsi process, if index too high set gpu 0
      Integer, Intent(In) :: index
      ! Local Variables
      Integer :: ierr, avail_devs, idx
      ierr = cudaGetDeviceCount(avail_devs)
      If(avail_devs.le.index) Then
        idx = 0
        Write(*,*) "**CUDA** ERROR: required device index too high, device 0 used instead"
      Else
        idx = index
      Endif
      ierr = cudaSetDevice(idx)
    End Subroutine bk_set_current_device

    Subroutine bk_initialize_gpus()
      Integer :: ierr
      If(do_DEV_cuda) Then
        If(init_blas_handle_status.eq.0) Then
          ierr = cublascreate(blas_handle)
          init_blas_handle_status = 1
        Endif
        If(init_solver_handle_status.eq.0) Then
          ierr = cusolverdncreate(solver_handle)
          init_solver_handle_status = 1
        Endif
        If(init_solver_params_status.eq.0) Then
          ierr = cusolverdncreateparams(solver_params)
          init_solver_params_status = 1
        Endif
      Endif
    End Subroutine bk_initialize_gpus

    Subroutine bk_mem_print(name)
    !! print in output memory usage of GPU currently associated with this process
    !! Only used for debug purposes
      Character(Len=12), Intent(In) :: name
      ! Local Variables
      Integer                       :: ierr, dev
      Integer(kind=cuda_count_kind) :: free, total
      Real(double)                  :: free_GB, total_GB
      ierr = cudaGetDevice(dev)
      ierr = cudaMemGetInfo(free,total)
      free_GB  = Real(free,Kind=double)/(1024.0_double**3)
      total_GB = Real(total,Kind=double)/(1024.0_double**3)
      Write(*,"(A,I1,A,F7.3,A,F7.3,A)") name//" =======CUDA=====> DEVICE ",dev," MEMORY:",total_GB-free_GB," / ",total_GB," GB"
    End Subroutine bk_mem_print

    Subroutine bk_get_memory_device(used_mem,total_mem)
    !! return used memory and total memory (MB) of the GPU currently
    !! associated with this process
      Real(double), Intent(Out) :: used_mem
      Real(double), Intent(Out) :: total_mem
      ! Local Variables
      Integer                       :: ierr
      Integer(Kind=cuda_count_kind) :: free_byte, total_byte
      ierr = cudaMemGetInfo(free_byte,total_byte)
      total_mem = Real(total_byte,Kind=double)/(1024.0_double**2)
      used_mem = total_mem - Real(free_byte,Kind=double)/(1024.0_double**2)
    End Subroutine bk_get_memory_device

!! --- Memory Managment ---
    Subroutine bk_real_alloc(a,n,m)
      Type(bk_real_matrix), Intent(  Out) :: a
      Integer             , Intent(In   ) :: n
      Integer             , Intent(In   ) :: m
      ! Local Variables
      Integer :: ierr
      If(do_HOST_lapack) Then
        Allocate(a%data(n,m), stat=ierr)
        Call error_check(ierr,"bk_real_alloc",-1)
        a%data = Huge(a%data)
      Endif
      If(do_DEV_cuda) Then
        Allocate(a%ddata(n,m), stat=ierr)
        Call error_check(ierr,"bk_real_alloc",1)
      Endif
    End Subroutine bk_real_alloc

    Subroutine bk_comp_alloc(a,n,m)
      Type(bk_comp_matrix), Intent(  Out) :: a
      Integer             , Intent(In   ) :: n
      Integer             , Intent(In   ) :: m
      ! Local Variables
      Integer :: ierr
      If(do_HOST_lapack) Then
        Allocate(a%data(n,m), stat=ierr)
        Call error_check(ierr,"bk_comp_alloc",-1)
        a%data = Huge(Real(a%data,Kind=double))
      Endif
      If(do_DEV_cuda) Then
        Allocate(a%ddata(n,m), stat=ierr)
        Call error_check(ierr,"bk_comp_alloc",1)
      Endif
    End Subroutine bk_comp_alloc

    Subroutine bk_real_dealloc(a)
      Type(bk_real_matrix), Intent(InOut) :: a
      ! Local Variables
      Integer :: ierr
      If(do_HOST_lapack) Then
        If(.not.Allocated(a%data)) Call error_check(-9,"bk_real_dealloc",-1)
        Deallocate(a%data, stat=ierr)
        Call error_check(ierr,"bk_real_dealloc",-2)
      Endif
      If(do_DEV_cuda) Then
        If(.not.Allocated(a%ddata)) Call error_check(-9,"bk_real_dealloc",1)
        Deallocate(a%ddata, stat=ierr)
        Call error_check(ierr,"bk_real_dealloc",2)
      Endif
    End Subroutine bk_real_dealloc

    Subroutine bk_comp_dealloc(a)
      Type(bk_comp_matrix), Intent(InOut) :: a
      ! Local Variables
      Integer :: ierr
      If(do_HOST_lapack) Then
        If(.not.Allocated(a%data)) Call error_check(-9,"bk_comp_dealloc",-1)
        Deallocate(a%data, stat=ierr)
        Call error_check(ierr,"bk_comp_dealloc",-2)
      Endif
      If(do_DEV_cuda) Then
        If(.not.Allocated(a%ddata)) Call error_check(-9,"bk_comp_dealloc",1)
        Deallocate(a%ddata, stat=ierr)
        Call error_check(ierr,"bk_comp_dealloc",2)
      Endif
    End Subroutine bk_comp_dealloc

    Subroutine bk_alloc_1d(vec,n)
      Type(bk_vector_1D), Intent(InOut) :: vec
      Integer           , Intent(In   ) :: n
      ! Local Variables
      Integer :: ierr
      If(do_HOST_lapack) Then
        If(Allocated(vec%data)) Call error_check(-9,"bk_alloc_1d",-1)
        Allocate(vec%data(n), stat=ierr)
        Call error_check(ierr,"bk_alloc_1d",-2)
      Endif
      If(do_DEV_cuda) Then
        If(Allocated(vec%ddata)) Call error_check(-9,"bk_alloc_1d",1)
        Allocate(vec%ddata(n), stat=ierr)
        Call error_check(ierr,"bk_alloc_1d",2)
      Endif
    End Subroutine bk_alloc_1d

    Subroutine bk_dealloc_1d(vec)
      Type(bk_vector_1D), Intent(InOut) :: vec
      ! Local Variables
      Integer :: ierr
      If(do_HOST_lapack) Then
        If(.not.Allocated(vec%data)) Call error_check(-9,"bk_dealloc_1d",-1)
        Deallocate(vec%data, stat=ierr)
        Call error_check(ierr,"bk_dealloc_1d",-2)
      Endif
      If(do_DEV_cuda) Then
        If(.not.Allocated(vec%ddata)) Call error_check(-9,"bk_dealloc_1d",1)
        Deallocate(vec%ddata, stat=ierr)
        Call error_check(ierr,"bk_dealloc_1d",2)
      Endif
    End Subroutine bk_dealloc_1d

!! --- Data Movement ---
    Subroutine bk_set_raw_real(a,raw_data)
      Type(bk_real_matrix),                 Intent(  Out) :: a
      Real(double)        , Dimension(:,:), Intent(In   ) :: raw_data
      If(do_HOST_lapack) Then
        a%data = raw_data
      Endif
      If(do_DEV_cuda) Then
        a%ddata = raw_data
      Endif
    End Subroutine bk_set_raw_real

    Subroutine bk_set_raw_comp(a,raw_data)
      Type(bk_comp_matrix),                 Intent(  Out) :: a
      Complex(double)     , Dimension(:,:), Intent(In   ) :: raw_data
      If(do_HOST_lapack) Then
        a%data = raw_data
      Endif
      If(do_DEV_cuda) Then
        a%ddata = raw_data
      Endif
    End Subroutine bk_set_raw_comp

    Subroutine bk_get_raw_real(a,raw_data)
      Type(bk_real_matrix),                              Intent(In   ) :: a
      Real(double)        , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      ! Local Variables
      Integer :: ldm
    !! ERROR :: CAN NOT RUN BOTH VERSION HERE
      If(do_HOST_lapack) Then
        ldm = Size(a%data, Dim=1)
        Allocate(raw_data(ldm,ldm))
        raw_data = a%data
      Endif
      If(do_DEV_cuda) Then
        ldm = Size(a%ddata, Dim=1)
        Allocate(raw_data(ldm,ldm))
        raw_data = a%ddata
      Endif
    End Subroutine bk_get_raw_real

    Subroutine bk_get_raw_comp(a,raw_data)
      Type(bk_comp_matrix),                              Intent(In   ) :: a
      Complex(double)     , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      ! Local Variables
      Integer :: ldm
    !! ERROR :: CAN NOT RUN BOTH VERSION HERE
      If(do_HOST_lapack) Then
        ldm = Size(a%data, Dim=1)
        Allocate(raw_data(ldm,ldm))
        raw_data = a%data
      Endif
      If(do_DEV_cuda) Then
        ldm = Size(a%ddata, Dim=1)
        Allocate(raw_data(ldm,ldm))
        raw_data = a%ddata
      Endif
    End Subroutine bk_get_raw_comp

    Subroutine bk_get_raw_real_vect(a,raw_vect,kt)
      Type(bk_real_matrix),                            Intent(In   ) :: a
      Real(double)        , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
      Integer             ,                            Intent(In   ) :: kt
      ! Local Variables
      Integer     ,                 Parameter   :: ONLY_REAL_KS = 1
      Integer     ,                 Parameter   :: BOTH_TYPES_KS = 2
      Integer                                   :: i, j
      Real(double), Dimension(:,:), Allocatable :: temp
      Integer                                   :: ldm
    !! ERROR :: CAN NOT RUN BOTH HERE
    !! WARNING :: GPU VERSION CAN BE IMPROVED
      If(do_HOST_lapack) Then
        ldm = Size(a%data, Dim=1)
        Select Case(kt)
        Case(ONLY_REAL_KS)
          Allocate(raw_vect(ldm**2))
        Case(BOTH_TYPES_KS)
          Allocate(raw_vect(2*ldm**2))
        End Select
        Do j = 1, ldm
          Do i = 1, ldm
            raw_vect((j-1)*ldm+i) = a%data(i,j)
          Enddo
        Enddo
      Endif
      If(do_DEV_cuda) Then
        ldm = Size(a%ddata, Dim=1)
        Allocate(temp(ldm,ldm))
        temp = a%ddata
        Select Case(kt)
        Case(ONLY_REAL_KS)
          Allocate(raw_vect(ldm**2))
        Case(BOTH_TYPES_KS)
          Allocate(raw_vect(2*ldm**2))
        End Select
        Do j = 1, ldm
          Do i = 1, ldm
            raw_vect((j-1)*ldm+i) = temp(i,j)
          Enddo
        Enddo
      Endif
    End Subroutine bk_get_raw_real_vect

    Subroutine bk_get_raw_comp_vect(a,raw_vect,kt)
      Type(bk_comp_matrix),               Intent(In   ) :: a
      Real(double)        , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
      Integer             ,               Intent(In   ) :: kt
      ! Local Variables
      Integer        ,                 Parameter   :: BOTH_TYPES_KS = 2
      Integer                                      :: i, j
      Complex(double), Dimension(:,:), Allocatable :: temp
      Integer                                      :: ldm
      If(kt.ne.BOTH_TYPES_KS) Call error_check(-9,"bk_get_raw_comp_vect")
      If(do_HOST_lapack) Then
        ldm = Size(a%data, Dim=1)
        Allocate(raw_vect(2*ldm**2))
        Do j = 1, ldm
          Do i = 1, ldm
            raw_vect( (j-1)*ldm*2 + (i-1)*2 + 1 ) = Real(a%data(i,j), Kind=double)
            raw_vect( (j-1)*ldm*2 + (i-1)*2 + 2 ) = Aimag(a%data(i,j))
          Enddo
        Enddo
      Endif
      If(do_DEV_cuda) Then
        ldm = Size(a%ddata, Dim=1)
        Allocate(temp(ldm,ldm))
        temp = A%ddata
        Allocate(raw_vect(2*ldm**2))
        Do j = 1, ldm
          Do i = 1, ldm
            raw_vect( (j-1)*ldm*2 + (i-1)*2 + 1 ) = Real(temp(i,j), Kind=double)
            raw_vect( (j-1)*ldm*2 + (i-1)*2 + 2 ) = Aimag(temp(i,j))
          Enddo
        Enddo
      Endif
    End Subroutine bk_get_raw_comp_vect

    Subroutine bk_set_raw_1d(vec,raw)
      Type(bk_vector_1D),               Intent(InOut) :: vec
      Real(double)      , Dimension(:), Intent(In   ) :: raw
      If(do_HOST_lapack) Then
        vec%data = raw
      Endif
      If(do_DEV_cuda) Then
        vec%ddata = raw
      Endif
    End Subroutine bk_set_raw_1d

    Subroutine bk_get_raw_1d(vec,raw)
      Type(bk_vector_1D),                            Intent(InOut) :: vec
      Real(double)      , Dimension(:), Allocatable, Intent(  Out) :: raw
      ! Local Variables
      Integer :: ldm
      If(do_HOST_lapack) Then
        ldm = Size(vec%data, Dim=1)
        Allocate(raw(ldm))
        raw = vec%data
      Endif
      If(do_DEV_cuda) Then
        ldm = Size(vec%ddata, Dim=1)
        Allocate(raw(ldm))
        raw = vec%ddata
      Endif
    End Subroutine bk_get_raw_1d

!! --- LA Operations ---
    Subroutine bk_real_cholesky(matrix)
      Type(bk_real_matrix), Intent(InOut) :: matrix
      ! Local Variables
      Integer                                         :: n
      Integer                                         :: ierr
      Integer                                         :: lwork
      Real(double), Dimension(:), Allocatable, Device :: work_d
      Integer     ,                            Device :: ierr_d
      If(do_HOST_lapack) Then
        n = Size(matrix%data, dim=1)
        Call dpotrf('u',n,matrix%data,n,ierr)
        Call error_check(ierr,"bk_real_cholesky",-1)
      Endif
      If(do_DEV_cuda) Then
        n = Size(matrix%ddata, dim=1)
        ierr = cusolverdndpotrf_buffersize(solver_handle,CUBLAS_FILL_MODE_UPPER,n,matrix%ddata,n,lwork)
        Call error_check(ierr,"bk_real_cholesky",1)
        Allocate(work_d(lwork),stat = ierr)
        Call error_check(ierr,"bk_real_cholesky",2)
        ierr = cusolverdndpotrf(solver_handle,CUBLAS_FILL_MODE_UPPER,n,matrix%ddata,n,work_d,lwork,ierr_d)
        Call error_check(ierr,"bk_real_cholesky",3)
        ierr = ierr_d
        Call error_check(ierr,"bk_real_cholesky",4)
      Endif
    End Subroutine bk_real_cholesky

    Subroutine bk_comp_cholesky(matrix)
      Type(bk_comp_matrix), Intent(InOut) :: matrix
      ! Local Variables
      Integer                                            :: n
      Integer                                            :: ierr
      Integer                                            :: lwork
      Complex(double), Dimension(:), Allocatable, Device :: work_d
      Integer        ,                            Device :: ierr_d
      If(do_HOST_lapack) Then
        n = Size(matrix%data, dim=1)
        Call zpotrf('u',n,matrix%data,n,ierr)
        Call error_check(ierr,"bk_comp_cholesky",-1)
      Endif
      If(do_DEV_cuda) Then
        n = Size(matrix%ddata, dim=1)
        ierr = cusolverdnzpotrf_buffersize(solver_handle,CUBLAS_FILL_MODE_UPPER,n,matrix%ddata,n,lwork)
        Call error_check(ierr,"bk_comp_cholesky",1)
        Allocate(work_d(lwork),stat = ierr)
        Call error_check(ierr,"bk_comp_cholesky",2)
        ierr = cusolverdnzpotrf(solver_handle,CUBLAS_FILL_MODE_UPPER,n,matrix%ddata,n,work_d,lwork,ierr_d)
        Call error_check(ierr,"bk_comp_cholesky",3)
        ierr = ierr_d
        Call error_check(ierr,"bk_comp_cholesky",4)
      Endif
    End Subroutine bk_comp_cholesky

    Subroutine bk_real_set_lower_to(matrix,val)
      Type(bk_real_matrix), Intent(InOut) :: matrix
      Real(double)        , Intent(In   ) :: val
      ! Local Variables
      Integer :: n, i, j
      Real(double), Dimension(:,:), Device, Allocatable :: temp

      If(do_HOST_lapack) Then
        n = Size(matrix%data, dim=1)
        do i = 1, n-1
           do j = i+1, n
              matrix%data(j,i) = val
           enddo
        enddo
      Endif
      If(do_DEV_cuda) Then
        n = Size(matrix%ddata, dim=1)
        Allocate(temp(n,n))
        temp = matrix%ddata
        !$cuf kernel do(2) <<< *, * >>>
        do i = 1, n-1
          do j = 2, n
            if (j.le.i) cycle
            temp(j,i) = val
          enddo
        enddo
        matrix%ddata = temp
      Endif
    End Subroutine bk_real_set_lower_to

    Subroutine bk_comp_set_lower_to(matrix,val)
      Type(bk_comp_matrix), Intent(InOut) :: matrix
      Complex(double)     , Intent(In   ) :: val
      ! Local Variables
      Integer :: n, i, j
      Complex(double), Dimension(:,:), Device, Allocatable :: temp
      If(do_HOST_lapack) Then
        n = Size(matrix%data, dim=1)
        do i = 1, n-1
           do j = i+1, n
              matrix%data(j,i) = val
           enddo
        enddo
      Endif
      If(do_DEV_cuda) Then
        n = Size(matrix%ddata, dim=1)
        Allocate(temp(n,n))
        temp = matrix%ddata
        !$cuf kernel do(2) <<< *, * >>>
        do i = 1, n-1
          do j = 2, n
            if (j.le.i) cycle
            temp(j,i) = val
          enddo
        enddo
        matrix%ddata = temp
      Endif
    End Subroutine bk_comp_set_lower_to

    Subroutine bk_real_invert(matrix)
      Type(bk_real_matrix), Intent(InOut) :: matrix
      ! Local Variables
      Integer                                            :: n
      Integer                                            :: ierr
      Integer(double)                                    :: dlwork, hlwork
      Integer(double)                                    :: one = 1_double
      Integer(1)     , Dimension(:), Allocatable         :: hwork
      Integer(1)     , Dimension(:), Allocatable, Device :: dwork_d
      Integer        ,                            Device :: istat_d
      If(do_HOST_lapack) Then
        n = Size(matrix%data, dim=1)
        Call dtrtri('u', 'n', n, matrix%data, n,ierr)
        Call error_check(ierr,"bk_real_invert",1)
      Endif
      If(do_DEV_cuda) Then
        n = Size(matrix%ddata, dim=1)
        ierr = cusolverdnxtrtri_buffersize(solver_handle,CUBLAS_FILL_MODE_UPPER, &
                                          & CUBLAS_DIAG_NON_UNIT,n, &
                                          & CUDADATATYPE(CUDA_R_64F),matrix%ddata, &
                                          & n,dlwork,hlwork)
        Call error_check(ierr,"bk_real_invert",2)
        dlwork = Max(one,dlwork) * 4
        Allocate(dwork_d(dlwork),stat=ierr)
        Call error_check(ierr,"bk_real_invert",3)
        hlwork = Max(one,hlwork) * 4
        Allocate(hwork(hlwork),stat=ierr)
        Call error_check(ierr,"bk_real_invert",4)
        ierr = cusolverdnxtrtri(solver_handle,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,n,&
                               & CUDADATATYPE(CUDA_R_64F),matrix%ddata,n,dwork_d,dlwork,  &
                               & hwork,hlwork,istat_d)
        Call error_check(ierr,"bk_real_invert",5)
        ierr = istat_d
        Call error_check(ierr,"bk_real_invert",6)
      Endif
    End Subroutine bk_real_invert

    Subroutine bk_comp_invert(matrix)
      Type(bk_comp_matrix), Intent(InOut) :: matrix
      ! Local Variables
      Integer :: n
      Integer :: ierr
      Integer(double)                                    :: dlwork, hlwork
      Integer(double)                                    :: one = 1_double
      Integer(1)     , Dimension(:), Allocatable         :: hwork
      Integer(1)     , Dimension(:), Allocatable, Device :: dwork_d
      Integer        ,                            Device :: istat_d
      If(do_HOST_lapack) Then
        n = Size(matrix%data, dim=1)
        Call ztrtri('u','n',n,matrix%data,n,ierr)
        Call error_check(ierr,"bk_comp_invert",1)
      Endif
      If(do_DEV_cuda) Then
        n = Size(matrix%ddata, dim=1)
        ierr = cusolverdnxtrtri_buffersize(solver_handle,CUBLAS_FILL_MODE_UPPER, &
                                          & CUBLAS_DIAG_NON_UNIT,n, &
                                          & cudadatatype(CUDA_C_64F),matrix%ddata, &
                                          & n,dlwork,hlwork)
        Call error_check(ierr,"bk_comp_invert",2)
        dlwork = Max(one,dlwork) * 4
        Allocate(dwork_d(dlwork),stat=ierr)
        Call error_check(ierr,"bk_comp_invert",3)
        hlwork = Max(one,hlwork) * 4
        Allocate(hwork(hlwork),stat=ierr)
        Call error_check(ierr,"bk_comp_invert",4)
        ierr = cusolverdnxtrtri(solver_handle,CUBLAS_FILL_MODE_UPPER,CUBLAS_DIAG_NON_UNIT,n, &
                               & cudadatatype(CUDA_C_64F),matrix%ddata,n,dwork_d,dlwork, &
                               & hwork,hlwork,istat_d)
        Call error_check(ierr,"bk_comp_invert",5)
        ierr = istat_d
        Call error_check(ierr,"bk_comp_invert",6)
      Endif
    End Subroutine bk_comp_invert

    Subroutine bk_real_mult(res,a,lopa,b,lopb)
    !! compurte res = alpha * op(a) * op(b) + beta * res
      Type(bk_real_matrix), Intent(InOut) :: res
      Type(bk_real_matrix), Intent(In   ) :: a
      Logical     ,         Intent(In   ) :: lopa
      Type(bk_real_matrix), Intent(In   ) :: b
      Logical     ,         Intent(In   ) :: lopb
      ! Local Variables
      Integer          :: ndim
      Character(Len=1) :: opa
      Character(Len=1) :: opb
      Integer          :: op_a
      Integer          :: op_b
      Real(double)     :: alpha = 1.0_double
      Real(double)     :: beta  = 0.0_double
      Integer          :: ierr
      If(do_HOST_lapack) Then
        ndim = max(Size(a%data,Dim=1),Size(b%data,Dim=1),Size(res%data,Dim=1))
        opa = Merge('t','n',lopa)
        opb = Merge('t','n',lopb)
        Call dgemm(opa,opb,ndim,ndim,ndim,alpha,a%data,ndim,b%data,ndim,beta,res%data,ndim)
      Endif
      If(do_DEV_cuda) Then
        ndim = max(Size(a%ddata,Dim=1),Size(b%ddata,Dim=1),Size(res%ddata,Dim=1))
        op_a = Merge(CUBLAS_OP_T,CUBLAS_OP_N,lopa)
        op_b = Merge(CUBLAS_OP_T,CUBLAS_OP_N,lopb)
        ierr = cublasdgemm(blas_handle,op_a,op_b,ndim,ndim,ndim,alpha,a%ddata,ndim,b%ddata,ndim,beta,res%ddata,ndim)
        Call error_check(ierr,"bk_real_mult")
      Endif
    End Subroutine bk_real_mult

    Subroutine bk_comp_mult(res,a,lopa,b,lopb)
    !! compurte res = alpha * op(a) * op(b)   +   beta * res
      Type(bk_comp_matrix), Intent(InOut) :: res
      Type(bk_comp_matrix), Intent(In   ) :: a
      Logical             , Intent(In   ) :: lopa
      Type(bk_comp_matrix), Intent(In   ) :: b
      Logical             , Intent(In   ) :: lopb
      ! Local Variables
      Integer          :: ndim
      Character(Len=1) :: opa
      Character(Len=1) :: opb
      Integer          :: op_a
      Integer          :: op_b
      Complex(double)  :: alpha = (1.0_double,0.0_double)
      Complex(double)  :: beta  = (0.0_double,0.0_double)
      Integer          :: ierr
      If(do_HOST_lapack) Then
        ndim = max(Size(a%data,Dim=1),Size(b%data,Dim=1),Size(res%data,Dim=1))
        opa = Merge('c','n',lopa)
        opb = Merge('c','n',lopb)
        Call zgemm(opa,opb,ndim,ndim,ndim,alpha,a%data,ndim,b%data,ndim,beta,res%data,ndim)
      Endif
      If(do_DEV_cuda) Then
        ndim = max(Size(a%ddata,Dim=1),Size(b%ddata,Dim=1),Size(res%ddata,Dim=1))
        op_a = Merge(CUBLAS_OP_C,CUBLAS_OP_N,lopa)
        op_b = Merge(CUBLAS_OP_C,CUBLAS_OP_N,lopb)
        ierr = cublaszgemm(blas_handle,op_a,op_b,ndim,ndim,ndim,alpha,a%ddata,ndim,b%ddata,ndim,beta,res%ddata,ndim)
        Call error_check(ierr,"bk_comp_mult")
      Endif
    End Subroutine bk_comp_mult

    Subroutine bk_real_sym_diag(matrix,evals)
    !! Compute all eigenvalues and eigenvectors of matrix
    !! matrix is lower triangular
      Type(bk_real_matrix),               Intent(InOut) :: matrix
      Type(bk_vector_1D), Intent(InOut) :: evals
      ! Local Variables
      Integer                                 :: ierr
      Integer                                 :: ldm
      Integer                                 :: lwork   !! dimension of work
      Integer                                 :: liwork  !! dimension of iwork
      Real(double), Dimension(:), Allocatable :: work
      Integer     , Dimension(:), Allocatable :: iwork
      Integer(double)                                    :: hlwork, dlwork
      Integer(1)     , Dimension(:), Allocatable         :: hwork
      Integer(1)     , Dimension(:), Allocatable, Device :: dwork_d
      Integer        ,                            Device :: istat_d
      If(do_HOST_lapack) Then
        ldm = Size(matrix%data,Dim=1)
        Allocate(work(1))
        Allocate(iwork(1))
        !! Inquiring call, just compute the optimal workspace
        Call dsyevd('V','L',ldm,matrix%data,ldm,evals%data,work,-1,iwork,-1,ierr)
        Call error_check(ierr,"bk_real_sym_diag",1)
        lwork  = work(1)
        liwork = iwork(1)
        Deallocate(work)
        Deallocate(iwork)
        Allocate(work(lwork))      !! Reallocate work
        Allocate(iwork(liwork))    !! Reallocate iwork
        !! Actual diagonalization
        Call dsyevd('V','L',ldm,matrix%data,ldm,evals%data,work,lwork,iwork,liwork,ierr)
        Call error_check(ierr,"bk_real_sym_diag",2)
      Endif
      If(do_DEV_cuda) Then
        ldm = Size(matrix%ddata,Dim=1)
        hlwork = -1_double
        dlwork = -1_double
        ierr = cusolverdnsetadvoptions(solver_params,0,CUSOLVER_ALG_0)
        Call error_check(ierr,"bk_real_sym_diag",3)
        ierr = cusolverdnxsyevd_buffersize(solver_handle,solver_params,CUSOLVER_EIG_MODE_VECTOR,    &
                                          & CUBLAS_FILL_MODE_LOWER,ldm, cudadatatype(CUDA_R_64F),   &
                                          & matrix%ddata,ldm, cudadatatype(CUDA_R_64F),evals%ddata, &
                                          & cudadatatype(CUDA_R_64F),dlwork,hlwork)
        Call error_check(ierr,"bk_real_sym_diag",4)
        Allocate(dwork_d(dlwork),stat=ierr)
        Call error_check(ierr,"bk_real_sym_diag",5)
        If(hlwork.le.0) hlwork = 1_double
        Allocate(hwork(hlwork),stat=ierr)
        Call error_check(ierr,"bk_real_sym_diag",6)
        ierr = cusolverdnxsyevd(solver_handle,solver_params,CUSOLVER_EIG_MODE_VECTOR,    &
                               & CUBLAS_FILL_MODE_LOWER,ldm, cudadatatype(CUDA_R_64F),   &
                               & matrix%ddata,ldm, cudadatatype(CUDA_R_64F),evals%ddata, &
                               & cudadatatype(CUDA_R_64F),dwork_d,dlwork,hwork,hlwork,istat_d)
        Call error_check(ierr,"bk_real_sym_diag",7)
        ierr = istat_d
        Call error_check(ierr,"bk_real_sym_diag",8)
        Deallocate(dwork_d)
        Deallocate(hwork)
      Endif
    End Subroutine bk_real_sym_diag

    Subroutine bk_comp_sym_diag(matrix,evals)
    !! Compute all eigenvalues and eigenvectors of matrix
    !! matrix is lower triangular
      Type(bk_comp_matrix),               Intent(InOut) :: matrix
      Type(bk_vector_1D), Intent(InOut) :: evals
      ! Local Variables
      Integer                                    :: ierr
      Integer                                    :: ldm
      Integer                                    :: lwork   !! dimension of work
      Integer                                    :: lrwork  !! dimension of rwork
      Integer                                    :: liwork  !! dimension of iwork
      Complex(double), Dimension(:), Allocatable :: work
      Real(double)   , Dimension(:), Allocatable :: rwork
      Integer        , Dimension(:), Allocatable :: iwork
      Integer(double)                                    :: hlwork, dlwork
      Integer(1)     , Dimension(:), Allocatable         :: hwork
      Integer(1)     , Dimension(:), Allocatable, Device :: dwork_d
      Integer        ,                            Device :: istat_d
      If(do_HOST_lapack) Then
        ldm = Size(matrix%data,Dim=1)
        Allocate(work(1))
        Allocate(rwork(1))
        Allocate(iwork(1))
        !! Inquiring call, just compute the optimal workspace
        Call zheevd('V','L',ldm,matrix%data,ldm,evals%data,work,-1,rwork,-1,iwork,-1,ierr)
        Call error_check(ierr,"bk_comp_sym_diag",1)
        lwork  = work(1)
        lrwork = rwork(1)
        liwork = iwork(1)
        Deallocate(work)
        Deallocate(rwork)
        Deallocate(iwork)
        Allocate(work(lwork))      !! Reallocate work
        Allocate(rwork(lrwork))    !! Reallocate rwork
        Allocate(iwork(liwork))    !! Reallocate iwork
        !! Actual diagonalization
        Call zheevd('V','L',ldm,matrix%data,ldm,evals%data,work,lwork,rwork,lrwork,iwork,liwork,ierr)
        Call error_check(ierr,"bk_comp_sym_diag",2)
      Endif
      If(do_DEV_cuda) Then
        ldm = Size(matrix%ddata,Dim=1)
        hlwork = -1_double
        dlwork = -1_double
        ierr = cusolverdnsetadvoptions(solver_params,0,CUSOLVER_ALG_0)
        Call error_check(ierr,"bk_comp_sym_diag",3)
        ierr = cusolverdnxsyevd_buffersize(solver_handle,solver_params,CUSOLVER_EIG_MODE_VECTOR,    &
                                          & CUBLAS_FILL_MODE_LOWER,ldm, cudadatatype(CUDA_C_64F),   &
                                          & matrix%ddata,ldm, cudadatatype(CUDA_R_64F),evals%ddata, &
                                          & cudadatatype(CUDA_C_64F),dlwork,hlwork)
        Call error_check(ierr,"bk_comp_sym_diag",4)
        Allocate(dwork_d(dlwork),stat=ierr)
        Call error_check(ierr,"bk_comp_sym_diag",5)
        If(hlwork.le.0) hlwork = 1_double
        Allocate(hwork(hlwork),stat=ierr)
        Call error_check(ierr,"bk_comp_sym_diag",6)
        ierr = cusolverdnxsyevd(solver_handle,solver_params,CUSOLVER_EIG_MODE_VECTOR,    &
                               & CUBLAS_FILL_MODE_LOWER,ldm, cudadatatype(CUDA_C_64F),   &
                               & matrix%ddata,ldm, cudadatatype(CUDA_R_64F),evals%ddata, &
                               & cudadatatype(CUDA_C_64F),dwork_d,dlwork,hwork,hlwork,istat_d)
        Call error_check(ierr,"bk_comp_sym_diag",7)
        ierr = istat_d
        Call error_check(ierr,"bk_comp_sym_diag",8)
        Deallocate(dwork_d)
        Deallocate(hwork)
      Endif
    End Subroutine bk_comp_sym_diag

    Subroutine bk_real_matrix_copy(a,b)
      Type(bk_real_matrix), Intent(InOut) :: a
      Type(bk_real_matrix), Intent(InOut) :: b
      ! Local Variables
      Integer :: ierr
      If(do_HOST_lapack) Then
        If(Allocated(a%data)) Deallocate(a%data)
        Call move_alloc(b%data, a%data)
      Endif
      If(do_DEV_cuda) Then
        Allocate(a%ddata(Size(b%ddata,Dim=1),Size(b%ddata,Dim=2)),stat=ierr)
        Call error_check(ierr,"bk_real_matrix_copy",1)
        a%ddata = b%ddata
        Deallocate(b%ddata,stat=ierr)
        Call error_check(ierr,"bk_real_matrix_copy",2)
      Endif
    End Subroutine bk_real_matrix_copy

    Subroutine bk_comp_matrix_copy(a,b)
      Type(bk_comp_matrix), Intent(InOut) :: a
      Type(bk_comp_matrix), Intent(InOut) :: b
      ! Local Variables
      Integer :: ierr
      If(do_HOST_lapack) Then
        If(Allocated(a%data)) Deallocate(a%data)
        Call move_alloc(b%data, a%data)
      Endif
      If(do_DEV_cuda) Then
        Allocate(a%ddata(Size(b%ddata,Dim=1),Size(b%ddata,Dim=2)),stat=ierr)
        Call error_check(ierr,"bk_comp_matrix_copy",1)
        a%ddata = b%ddata
        Deallocate(b%ddata,stat=ierr)
        Call error_check(ierr,"bk_comp_matrix_copy",2)
      Endif
    End Subroutine bk_comp_matrix_copy

!! --- Matrix Related ---
    Function bk_real_allocated(matrix) Result(res)
      Logical                          :: res
      Type(bk_real_matrix), Intent(In) :: matrix
      If(do_HOST_lapack) Then
        res = Allocated(matrix%data)
      Endif
      If(do_DEV_cuda) Then
        res = Allocated(matrix%ddata)
      Endif
    End Function bk_real_allocated

    Function bk_comp_allocated(matrix) Result(res)
      Logical                          :: res
      Type(bk_comp_matrix), Intent(In) :: matrix
      If(do_HOST_lapack) Then
        res = Allocated(matrix%data)
      Endif
      If(do_DEV_cuda) Then
        res = Allocated(matrix%ddata)
      Endif
    End Function bk_comp_allocated

    Subroutine bk_real_shift_diag(matrix,shift,nocc)
      Type(bk_real_matrix), Intent(InOut) :: matrix
      Real(double)        , Intent(In   ) :: shift
      Integer             , Intent(In   ) :: nocc
      ! Local Variables
      Integer :: ldm
      Integer :: m
      Integer :: i, j
      Real(double), Dimension(:,:), Device, Allocatable :: temp
      If(do_HOST_lapack) Then
        ldm = Size(matrix%data,Dim=1)
        Loop_on_diag: Do m = 1, nocc
          matrix%data(m,m) = matrix%data(m,m) - shift
        Enddo Loop_on_diag
      Endif
      If(do_DEV_cuda) Then
        ldm = Size(matrix%ddata,Dim=1)
        Allocate(temp(ldm,ldm))
        temp = matrix%ddata
        !$cuf kernel do(1) <<< *, * >>>
        Do i = 1, nocc
          temp(i,i) = temp(i,i) - shift
        Enddo
        matrix%ddata= temp
      Endif
    End Subroutine bk_real_shift_diag

    Subroutine bk_comp_shift_diag(matrix,shift,nocc)
    !! Just shift the real component on the diagonal
      Type(bk_comp_matrix), Intent(InOut) :: matrix
      Real(double)        , Intent(In   ) :: shift
      Integer             , Intent(In   ) :: nocc
      ! Local Variables
      Integer :: ldm
      Integer :: m
      Integer :: i, j
      Complex(double), Dimension(:,:), Device, Allocatable :: temp
      If(do_HOST_lapack) Then
        ldm = Size(matrix%data,Dim=1)
        Loop_on_diag: Do m = 1, nocc
          matrix%data(m,m) = matrix%data(m,m) - Cmplx(shift,0.0_double,Kind=double)
        Enddo Loop_on_diag
      Endif
      If(do_DEV_cuda) Then
        ldm = Size(matrix%ddata,Dim=1)
        Allocate(temp(ldm,ldm))
        temp = matrix%ddata
        !$cuf kernel do(1) <<< *, * >>>
        Do i = 1, nocc
          temp(i,i) = temp(i,i) - Cmplx(shift,0.0_double,Kind=double)
        Enddo
        matrix%ddata= temp
      Endif
    End Subroutine bk_comp_shift_diag

    Subroutine bk_compute_real_tester_sub(matrix,nocc,tstr)
      Type(bk_real_matrix), Intent(In   ) :: matrix
      Integer             , Intent(In   ) :: nocc
      Real(double)        , Intent(  Out) :: tstr
      ! Local Variables
      Integer                                 :: ldm
      Integer                                 :: m
      Integer                                 :: n
      Real(double), Dimension(:), Allocatable :: diagonal
      If(do_HOST_lapack) Then
        ldm = Size(matrix%data,Dim=1)
        Allocate(diagonal(ldm))
        Do n = 1, ldm
          diagonal(n) = matrix%data(n,n)
        Enddo
        tstr = 0.0_double
        Loop_on_cols: Do m = 1, nocc
          Loop_on_rows: Do n = nocc+1, ldm
            Associate( ediff => Abs( diagonal(m)-diagonal(n) ) + 0.004_double )
              tstr = Max( matrix%data(n,m)**2/ediff , tstr )
            End Associate
          Enddo Loop_on_rows
        Enddo Loop_on_cols
      Endif
      If(do_DEV_cuda) Then
        tstr = -66.0_double
      Endif
      !!$cuf kernel do(2) <<< *, * >>>
      !do i = 1, ldm-1
      !  do j = 2, ldm
      !    if (j.le.i) cycle
      !    temp(j,i) = 0.0_double
      !  enddo
      !enddo
    End Subroutine bk_compute_real_tester_sub

    Subroutine bk_compute_complex_tester_sub(matrix,nocc,tstr)
    !! Just computing tester value on the real part of complex matrices
      Type(bk_comp_matrix), Intent(In   ) :: matrix
      Integer             , Intent(In   ) :: nocc
      Real(double)        , Intent(  Out) :: tstr
      ! Local Variables
      Integer                                 :: ldm
      Integer                                 :: m
      Integer                                 :: n
      Real(double), Dimension(:), Allocatable :: diagonal
      If(do_HOST_lapack) Then
        ldm = Size(matrix%data,Dim=1)
        Allocate(diagonal(ldm))
        Do n = 1, ldm
          diagonal(n) = Real(matrix%data(n,n))
        Enddo
        tstr = 0.0_double
        Loop_on_cols: Do m = 1, nocc
          Loop_on_rows: Do n = nocc+1, ldm
            Associate( ediff => Abs( diagonal(m)-diagonal(n) ) + 0.004_double )
              tstr = Max( Real(matrix%data(n,m))**2/ediff , tstr )
            End Associate
          Enddo Loop_on_rows
        Enddo Loop_on_cols
      Endif
      If(do_DEV_cuda) Then
        tstr = -66.0_double
      Endif
    End Subroutine bk_compute_complex_tester_sub

    ! Private
    Subroutine error_check(code,name,step)
      Integer         , Intent(In)           :: code
      Character(Len=*), Intent(In)           :: name
      Integer         , Intent(In), Optional :: step
      If(code.ne.0) Then
        If(Present(step)) Then
          Write(*,"(A,I0,3X,'(',I0,')')") "**CUDA** ERROR in "//name//" step: ",step,code
        Else
          Write(*,"(A,3X,'(',I0,')')") "**CUDA** ERROR in "//name,code
        Endif
      Endif
    End Subroutine error_check

End Module backend_module
