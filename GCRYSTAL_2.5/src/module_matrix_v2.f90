Module module_matrix_v2
!===========================================================
!    @author: Giacomo Ambrogio       date: Feb 2024
!
!    DASH - Device Accelerated Solutions for HPC
!
!    Matrix module class
!
!    This module handles the matrices, it is the innermost
!    that do not depends on device languages
!
!    It uses abstract to deal with real/complex datatypes
!===========================================================
Use base_numbers
Use backend_module

  Implicit None
  Public

  Type, Public :: ks_vector_1D
  !! ============================================================
  !! Vector type (1D) used to deal with eigenvalues, only real
  !! ============================================================
    Type(bk_vector_1D), Private :: data
    Contains
    Procedure, Public :: create          => ks_vector_1D_create
    Procedure, Public :: destroy         => ks_vector_1D_destroy
    Procedure, Public :: get_raw         => ks_vector_1D_get_raw
    Procedure, Public :: set_raw         => ks_vector_1D_set_raw
    Procedure, Public :: shift           => ks_vector_1D_shift
    Procedure, Public :: sort            => ks_vector_1D_sort
    Procedure, Public :: get_bands_range => ks_vector_1D_get_bands_range
  End Type ks_vector_1D

  Type, Abstract, Public :: ks_matrix
  !! ============================================================
  !! Abstract type to deal with real and complex ks points
  !! This is integrated with the Extends(ks_matrix) types below
  !! Each of the extends types provides an implementation of the
  !! deferred procedures listed here
  !! ============================================================
  Logical, Private :: daggered = .False.
    Contains
    ! Public methods and interfaces
    Procedure(create)         , Deferred, Public :: create
    Procedure(free)           , Deferred, Public :: free
    Generic                   ,           Public :: set_raw        => set_raw_real     , set_raw_complex
    Generic                   ,           Public :: get_raw        => get_raw_real     , get_raw_complex     !, get_raw_vect
    Generic                   ,           Public :: host_set_raw   => host_set_raw_real, host_set_raw_complex
    Generic                   ,           Public :: host_get_raw   => host_get_raw_real, host_get_raw_complex
    Procedure(myallocated)    , Deferred, Public :: allocated
    Procedure(diag)           , Deferred, Public :: diag
    Procedure(copy)           , Deferred, Public :: copy
    Procedure(compute_tester) , Deferred, Public :: compute_tester
    Procedure(shift_diagonal) , Deferred, Public :: shift_diagonal
    Procedure                 ,           Public :: dagger         => ks_matrix_set_dagger
    Procedure(binary_op)      , Deferred, Public :: multiply
    Procedure(unary_op)       , Deferred, Public :: cholesky
    Procedure(unary_op)       , Deferred, Public :: clean_lower
    Procedure(unary_op)       , Deferred, Public :: invert
    Procedure(create)         , Deferred, Public :: host_create
    Procedure(free)           , Deferred, Public :: host_free
    Procedure(copy)           , Deferred, Public :: host_copy
    Procedure(unary_op)       , Deferred, Public :: host_to_dev
    Procedure(unary_op)       , Deferred, Public :: dev_to_host
    Procedure(dens)           , Deferred, Public :: compute_density
    ! Private methods
    Procedure(set_raw_real)   , Deferred, Private :: set_raw_real
    Procedure(set_raw_complex), Deferred, Private :: set_raw_complex
    Procedure(get_raw_real)   , Deferred, Private :: get_raw_real
    Procedure(get_raw_complex), Deferred, Private :: get_raw_complex
    Procedure(set_raw_real)   , Deferred, Private :: host_set_raw_real
    Procedure(set_raw_complex), Deferred, Private :: host_set_raw_complex
    Procedure(get_raw_real)   , Deferred, Private :: host_get_raw_real
    Procedure(get_raw_complex), Deferred, Private :: host_get_raw_complex
  End type ks_matrix

  Type, Extends(ks_matrix), Public :: real_ks_matrix
  !! ============================================================
  !! Extension of ks_matrix for real ks points
  !! It provides type-specific implementations:
  !!   basically manage type-specific calls to backend procedures
  !! ============================================================
    Type(bk_real_matrix), Private :: data
    Contains
    ! Public methods
    Procedure, Public :: create          => real_matrix_create
    Procedure, Public :: free            => real_matrix_free
    Procedure, Public :: allocated       => real_matrix_allocated
    Procedure, Public :: diag            => real_matrix_diag
    Procedure, Public :: copy            => real_matrix_copy
    Procedure, Public :: compute_tester  => real_matrix_compute_tester
    Procedure, Public :: shift_diagonal  => real_matrix_shift_diagonal
    Procedure, Public :: multiply        => real_matrix_multiply
    Procedure, Public :: cholesky        => real_matrix_cholesky
    Procedure, Public :: clean_lower     => real_matrix_clean_lower
    Procedure, Public :: invert          => real_matrix_invert
    Procedure, Public :: host_create     => real_matrix_host_create
    Procedure, Public :: host_free       => real_matrix_host_free
    Procedure, Public :: host_copy       => real_matrix_host_copy
    Procedure, Public :: host_to_dev     => real_matrix_host_to_dev
    Procedure, Public :: dev_to_host     => real_matrix_dev_to_host
    Procedure, Public :: compute_density => real_matrix_compute_density
    ! Private implementations
    Procedure, Private :: set_raw_real         => real_matrix_set_raw_real
    Procedure, Private :: set_raw_complex      => real_matrix_set_raw_complex
    Procedure, Private :: get_raw_real         => real_matrix_get_raw_real
    Procedure, Private :: get_raw_complex      => real_matrix_get_raw_complex
    Procedure, Private :: host_set_raw_real    => real_matrix_host_set_raw_real
    Procedure, Private :: host_set_raw_complex => real_matrix_host_set_raw_complex
    Procedure, Private :: host_get_raw_real    => real_matrix_host_get_raw_real
    Procedure, Private :: host_get_raw_complex => real_matrix_host_get_raw_complex
  End type real_ks_matrix
!
  Type, Extends(ks_matrix), Public :: complex_ks_matrix
  !! ============================================================
  !! Extension of ks_matrix for complex ks points
  !! It provides type-specific implementations:
  !!   basically manage type-specific calls to backend procedures
  !! ============================================================
    Type(bk_comp_matrix), Private :: data
    Contains
    ! Public methods
    Procedure, Public :: create          => complex_matrix_create
    Procedure, Public :: free            => complex_matrix_free
    Procedure, Public :: allocated       => complex_matrix_allocated
    Procedure, Public :: diag            => complex_matrix_diag
    Procedure, Public :: copy            => complex_matrix_copy
    Procedure, Public :: compute_tester  => complex_matrix_compute_tester
    Procedure, Public :: shift_diagonal  => complex_matrix_shift_diagonal
    Procedure, Public :: multiply        => complex_matrix_multiply
    Procedure, Public :: cholesky        => complex_matrix_cholesky
    Procedure, Public :: clean_lower     => complex_matrix_clean_lower
    Procedure, Public :: invert          => complex_matrix_invert
    Procedure, Public :: host_create     => complex_matrix_host_create
    Procedure, Public :: host_free       => complex_matrix_host_free
    Procedure, Public :: host_copy       => complex_matrix_host_copy
    Procedure, Public :: host_to_dev     => complex_matrix_host_to_dev
    Procedure, Public :: dev_to_host     => complex_matrix_dev_to_host
    Procedure, Public :: compute_density => complex_matrix_compute_density
    ! Private implementations
    Procedure, Private :: set_raw_real         => complex_matrix_set_raw_real
    Procedure, Private :: set_raw_complex      => complex_matrix_set_raw_complex
    Procedure, Private :: get_raw_real         => complex_matrix_get_raw_real
    Procedure, Private :: get_raw_complex      => complex_matrix_get_raw_complex
    Procedure, Private :: host_set_raw_real    => complex_matrix_host_set_raw_real
    Procedure, Private :: host_set_raw_complex => complex_matrix_host_set_raw_complex
    Procedure, Private :: host_get_raw_real    => complex_matrix_host_get_raw_real
    Procedure, Private :: host_get_raw_complex => complex_matrix_host_get_raw_complex
  End type complex_ks_matrix

  Abstract Interface
  !! ============================================================
  !! Big interface block necerrasy for the ks_matrix abstract
  !! derived type
  !! Here we need an explicit interface for all the deferred
  !! procedure listed in ks_matrix (and implemehted in the exteds
  !! real_ks_matrix and complex_ks_matrix types)
  !! ============================================================
    Function myallocated(M) Result(res)
      Import :: ks_matrix
      Implicit None
      Logical :: res
      Class(ks_matrix), Intent(In) :: M
    End Function myallocated
    Subroutine create(M,n)
      Import :: ks_matrix
      Implicit None
      Class(ks_matrix), Intent(InOut) :: M
      Integer         , Intent(In   ) :: n
    End Subroutine create
    Subroutine free(M)
      Import :: ks_matrix
      Implicit None
      Class(ks_matrix), Intent(InOut) :: M
    End Subroutine free
    Function compute_tester(M,nocc) Result(tstr)
      Import :: fp64
      Import :: ks_matrix
      Real(fp64) :: tstr
      Class(ks_matrix), Intent(In) :: M
      Integer         , Intent(In) :: nocc
    End Function compute_tester
    Subroutine shift_diagonal(M,shift,nocc)
      Import :: fp64
      Import :: ks_matrix
      Class(ks_matrix), Intent(InOut) :: M
      Real(fp64)    , Intent(In) :: shift
      Integer         , Intent(In) :: nocc
    End Subroutine shift_diagonal
    Subroutine get_raw_real(M,raw_data)
      Import :: fp64
      Import :: ks_matrix
      Implicit None
      Class(ks_matrix),                              Intent(In   ) :: M
      Real(fp64)    , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
    End Subroutine get_raw_real
    !Subroutine get_raw_vect(M,raw_vect,kt)
    !  Import :: fp64
    !  Import :: ks_matrix
    !  Implicit None
    !  Class(ks_matrix),                            Intent(In   ) :: M
    !  Real(fp64)    , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
    !  Integer         ,                            Intent(In   ) :: kt
    !End Subroutine get_raw_vect
    Subroutine get_raw_complex(M,raw_data)
      Import :: fp64
      Import :: ks_matrix
      Implicit None
      Class(ks_matrix),                              Intent(In   ) :: M
      Complex(fp64) , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
    End Subroutine get_raw_complex
    Subroutine set_raw_real(M,raw_data)
      Import :: fp64
      Import :: ks_matrix
      Implicit None
      Class(ks_matrix),                 Intent(InOut) :: M
      Real(fp64)    , Dimension(:,:), Intent(In   ) :: raw_data
    End Subroutine set_raw_real
    Subroutine set_raw_complex(M,raw_data)
      Import :: fp64
      Import :: ks_matrix
      Implicit None
      Class(ks_matrix),                 Intent(InOut) :: M
      Complex(fp64) , Dimension(:,:), Intent(In   ) :: raw_data
    End Subroutine set_raw_complex
    Subroutine diag(M,evals)
      Import :: ks_matrix
      Import :: ks_vector_1D
      Class(ks_matrix)  , Intent(InOut) :: M
      Type(ks_vector_1D), Intent(InOut) :: evals
    End Subroutine diag
    Subroutine copy(M,N)
      Import :: ks_matrix
      Class(ks_matrix), Intent(InOut) :: M
      Class(ks_matrix), Intent(InOut) :: N
    End Subroutine copy
    !#######################################################
    Subroutine dens(M,nocc,pg,weights,coords)
      Import :: ks_matrix
      Import :: fp64
      Implicit None
      Class(ks_matrix),               Intent(InOut) :: M
      Integer         ,               Intent(In   ) :: nocc
      Real(fp64)    , Dimension(:), Intent(InOut) :: pg
      Real(fp64)    , Dimension(:), Intent(In   ) :: weights
      Integer         , Dimension(3), Intent(In   ) :: coords
    End Subroutine dens
    !#######################################################
    Subroutine unary_op(M)
    !! A unary operation on a base class object
      Import :: ks_matrix
      Implicit None
      Class(ks_matrix), Intent(InOut)  :: M
    End Subroutine unary_op
    Subroutine binary_op(M,N,O)
      !! A binary operation on a base class object
      Import :: ks_matrix
      Implicit None
      Class(ks_matrix), Intent(InOut) :: M
      Class(ks_matrix), Intent(In   ) :: N
      Class(ks_matrix), Intent(In   ) :: O
    End Subroutine binary_op
  End Interface

  Contains
    !! ============================================================
    !!        1D VECTOR PROCEDURES
    !! ============================================================
    Subroutine ks_vector_1D_create(V1D,ldp)
      Class(ks_vector_1D), Intent(InOut) :: V1D
      Integer           , Intent(In   ) :: ldp
      Call bk_alloc_1d(V1D%data,ldp)
    End Subroutine ks_vector_1D_create

    Subroutine ks_vector_1D_destroy(V1D)
      Class(ks_vector_1D), Intent(InOut) :: V1D
      Call bk_dealloc_1d(V1D%data)
    End Subroutine ks_vector_1D_destroy

    Subroutine ks_vector_1D_get_raw(V1D,raw_data)
      Class(ks_vector_1D),                      Intent(InOut) :: V1D
      Real(fp64), Dimension(:), Allocatable, Intent(  Out) :: raw_data
      Call bk_get_raw_1d(V1D%data,raw_data)
    End Subroutine ks_vector_1D_get_raw

    Subroutine ks_vector_1D_set_raw(V1D,raw_data)
    Use dash_utils, Only: dash_error
      Class(ks_vector_1D),         Intent(InOut) :: V1D
      Real(fp64), Dimension(:), Intent(In   ) :: raw_data
      Call dash_error(0,'ks_vector_1D_set_raw','NOT IMPLEMENTED')
      !Call bk_set_raw_1d(P1D%vector,raw_data)
    End Subroutine ks_vector_1D_set_raw

    Subroutine ks_vector_1D_shift(V1D,occ,shift)
      Class(ks_vector_1D), Intent(InOut) :: V1D
      Integer            , Intent(In   ) :: occ
      Real(fp64)       , Intent(In   ) :: shift
      Call bk_shift_1d(V1D%data,occ,shift)
    End Subroutine ks_vector_1D_shift

    Subroutine ks_vector_1D_sort(V1D,M)
      Class(ks_vector_1D), Intent(InOut)           :: V1D
      Class(ks_matrix)   , Intent(InOut), Optional :: M
      If(Present(M)) Then
        Select Type(M)
        Type Is(real_ks_matrix)
          Call bk_sort_eval_evecs(V1D%data,M%data)
        Type Is(complex_ks_matrix)
          Call bk_sort_eval_evecs(V1D%data,M%data)
        End Select
      Else
        Call bk_sort_eval_evecs(V1D%data)
      Endif
    End Subroutine ks_vector_1D_sort

    Subroutine ks_vector_1D_get_bands_range(V1D,this_k,this_s)
      Class(ks_vector_1D),               Intent(InOut) :: V1D
      Integer            ,               Intent(In   ) :: this_k
      Integer            ,               Intent(In   ) :: this_s
      !Integer            , Dimension(:), Intent(InOut) :: bmaxk
      !Integer            , Dimension(:), Intent(InOut) :: bmink
      !Real(fp64)       , Dimension(:), Intent(InOut) :: bmaxv
      !Real(fp64)       , Dimension(:), Intent(InOut) :: bminv
      !Call bk_get_bands_range_1d(V1D%data,this_k,this_s)
    End Subroutine ks_vector_1D_get_bands_range

    !! ============================================================
    !!        1D VECTOR PROCEDURES - END
    !! ============================================================

    Subroutine ks_matrix_set_dagger(M,op)
    !! ============================================================
    !! Set dagger falg to the ks_matrix M, can be improved...
    !! ============================================================
      Class(ks_matrix), Intent(InOut) :: M
      Integer         , Intent(In   ) :: op
      Select Case(op)
      Case(OP_NONE)
        M%daggered = .False.
      Case(OP_TRANS, OP_DAGGER)
        M%daggered = .True.
      End Select
    End Subroutine ks_matrix_set_dagger

    Function real_matrix_allocated(M) Result(res)
    !! ============================================================
    !! Allocated function implementation for the ks_matrix
    !! real version
    !! ============================================================
    Logical :: res
    Class(real_ks_matrix), Intent(In) :: M
      res = bk_real_allocated(M%data)
    End Function real_matrix_allocated

    Function complex_matrix_allocated(M) Result(res)
    !! ============================================================
    !! Allocated function implementation for the ks_matrix
    !! complex version
    !! ============================================================
    Logical :: res
    Class(complex_ks_matrix), Intent(In) :: M
      res = bk_comp_allocated(M%data)
    End Function complex_matrix_allocated





















    !! --- Store matrix in memory ---
    Subroutine real_matrix_create(M,n)
    !! Create the data for a real n*n matrix
      Class(real_ks_matrix ), Intent(InOut) :: M
      Integer               , Intent(In   ) :: n
      ! Local Variables
      Integer :: ldm
      ldm = max(1,n)
      Call bk_real_alloc(M%data,ldm,ldm)
    End Subroutine real_matrix_create
    Subroutine complex_matrix_create(M,n)
    !! Create the data for a complex n*n matrix
      Class(complex_ks_matrix), Intent(InOut) :: M
      Integer                 , Intent(In   ) :: n
      !Local Variables
      Integer :: ldm
      ldm = max(1,n)
      Call bk_comp_alloc(M%data,ldm,ldm)
    End Subroutine complex_matrix_create

    !! --- Free matrix memory ---
    Subroutine real_matrix_free(M)
    !! Create the data for a real n*n matrix
      Class(real_ks_matrix ), Intent(InOut) :: M
      Call bk_real_dealloc(M%data)
    End Subroutine real_matrix_free
    Subroutine complex_matrix_free(M)
    !! Create the data for a complex n*n matrix
      Class(complex_ks_matrix), Intent(InOut) :: M
      Call bk_comp_dealloc(M%data)
    End Subroutine complex_matrix_free

    !! --- Set raw data ---
    Subroutine real_matrix_set_raw_real(M,raw_data)
    !! Set the raw data for M
      Class(real_ks_matrix),                 Intent(InOut) :: M
      Real(fp64)         , Dimension(:,:), Intent(In   ) :: raw_data
      Call bk_set_raw_real(M%data,raw_data)
    End Subroutine real_matrix_set_raw_real
    Subroutine real_matrix_set_raw_complex(M,raw_data)
    !! Set the raw data for M
    Use dash_utils
      Class(real_ks_matrix),               Intent(InOut) :: M
      Complex(fp64)      , Dimension(:,:), Intent(In   ) :: raw_data
      Call dash_error(0,'real_matrix_set_raw_complex','invalid')
    End Subroutine real_matrix_set_raw_complex
    Subroutine complex_matrix_set_raw_real(M,raw_data)
    !! set the raw data for M
    Use dash_utils
      Class(complex_ks_matrix)                , Intent(InOut) :: M
      Real(fp64)            , Dimension(:,:), Intent(In   ) :: raw_data
      Call dash_error(0,'complex_matrix_set_raw_real','invalid assign')
      ! It is possible to assign a real matrix to a complex one,
      ! but still it is not the cleanest way to do it
      !M%data = raw_data
    End Subroutine complex_matrix_set_raw_real
    Subroutine complex_matrix_set_raw_complex(M,raw_data)
    !! Set the raw data for M
      Class(complex_ks_matrix)                , Intent(InOut) :: M
      Complex(fp64)         , Dimension(:,:), Intent(In   ) :: raw_data
      Call bk_set_raw_comp(M%data,raw_data)
    End Subroutine complex_matrix_set_raw_complex

    !! --- Get the raw data ---
    Subroutine real_matrix_get_raw_real(M,raw_data)
    !! Get the raw data for M
      Class(real_ks_matrix),                              Intent(In   ) :: M
      Real(fp64)         , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call bk_get_raw_real(M%data,raw_data)
    End Subroutine real_matrix_get_raw_real
    !Subroutine real_matrix_get_raw_vector(M,raw_vect,kt)
    !!! get raw data from M, put in CRYSTAL like 1D vector, std data ordering
    !  Class(real_ks_matrix),                            Intent(In   ) :: M
    !  Real(fp64)         , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
    !  Integer              ,                            Intent(In   ) :: kt
    !  Call bk_get_raw_real_vect(M%data,raw_vect,kt)
    !End Subroutine real_matrix_get_raw_vector
    Subroutine real_matrix_get_raw_complex(M,raw_data)
    !! Get the raw data for A
    Use dash_utils
      Class(real_ks_matrix),                                Intent(In   ) :: M
      Complex(fp64)        , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call dash_error(0,'real_matrix_get_raw_complex','invalid')
    End Subroutine real_matrix_get_raw_complex
    Subroutine complex_matrix_get_raw_real(M,raw_data)
    !! Get the raw data for M
    Use dash_utils
      Class(complex_ks_matrix)                ,              Intent(In   ) :: M
      Real(fp64)            , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call dash_error(0,'complex_matrix_get_raw_complex','invalid')
    End Subroutine complex_matrix_get_raw_real
    !Subroutine complex_matrix_get_raw_vector(M,raw_vect,kt)
    !!! Get the raw data for M
    !Use dash_utils
    !  Class(complex_ks_matrix),                            Intent(In   ) :: M
    !  Real(fp64)            , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
    !  Integer                 ,                            Intent(In   ) :: kt
    !  ! Local Variables
    !  Integer, Parameter :: ONLY_REAL_KS = 1
    !  Integer, Parameter :: BOTH_TYPES_KS = 2
    !  Integer            :: i, j
    !  Call bk_get_raw_comp_vect(M%data,raw_vect,kt)
    !End Subroutine complex_matrix_get_raw_vector
    Subroutine complex_matrix_get_raw_complex(M,raw_data)
    !! Get the raw data for M
      Class (complex_ks_matrix)                             , Intent(In   ) :: M
      Complex(fp64)          , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call bk_get_raw_comp(M%data,raw_data)
    End Subroutine complex_matrix_get_raw_complex

    Function real_matrix_compute_tester(M,nocc) Result(tstr)
      Real(fp64)                      :: tstr
      Class(real_ks_matrix), Intent(In) :: M
      Integer              , Intent(In) :: nocc

      Call bk_compute_real_tester_sub(M%data,nocc,tstr)
    End Function real_matrix_compute_tester

    Function complex_matrix_compute_tester(M,nocc) Result(tstr)
      Real(fp64)                         :: tstr
      Class(complex_ks_matrix), Intent(In) :: M
      Integer                 , Intent(In) :: nocc

      Call bk_compute_complex_tester_sub(M%data,nocc,tstr)
    End Function complex_matrix_compute_tester

    Subroutine real_matrix_shift_diagonal(M,shift,nocc)
      Class(real_ks_matrix), Intent(InOut) :: M
      Real(fp64)         , Intent(In   ) :: shift
      Integer              , Intent(In   ) :: nocc
      Call bk_real_shift_diag(M%data,shift,nocc)
    End Subroutine real_matrix_shift_diagonal

    Subroutine complex_matrix_shift_diagonal(M,shift,nocc)
      Class(complex_ks_matrix), Intent(InOut) :: M
      Real(fp64)            , Intent(In   ) :: shift
      Integer                 , Intent(In   ) :: nocc
      Call bk_comp_shift_diag(M%data,shift,nocc)
    End Subroutine complex_matrix_shift_diagonal

    Subroutine real_matrix_multiply(M,N,O)
    !! Multiply M and N matreces, considering their transposition status
    Use dash_utils
      Class(real_ks_matrix), Intent(InOut) :: M
      Class(ks_matrix)     , Intent(InOut) :: N
      Class(ks_matrix)     , Intent(InOut) :: O
      Select Type(N)
      Class Is(real_ks_matrix)
        Select Type(O)
        Class Is(real_ks_matrix)
          !! compute: M = N * O, daggered is to transpose the relative matrix
        Call bk_real_mult(M%data,N%data,N%daggered,O%data,O%daggered)
        Class Default
          Call dash_error(0,'real_matrix_multiply','invalid type combination')
        End Select
      Class Default
        Call dash_error(0,'real_matrix_multiply','invalid type combination')
      End Select
      M%daggered = .False.
      N%daggered = .False.
      O%daggered = .False.
    End Subroutine real_matrix_multiply
    Subroutine complex_matrix_multiply(M,N,O)
    !! Multiply M and N matreces, considering their transposition status
    Use dash_utils
      Class(complex_ks_matrix), Intent(InOut) :: M
      Class(ks_matrix)        , Intent(InOut) :: N
      Class(ks_matrix)        , Intent(InOut) :: O
      Select Type(N)
      Class Is(complex_ks_matrix)
        Select Type(O)
        Class Is(complex_ks_matrix)
          !! compute: M = N * O, daggered is to conj-transpose the relative matrix
          Call bk_comp_mult(M%data,N%data,N%daggered,O%data,O%daggered)
        Class Default
          Call dash_error(0,'complex_matrix_multiply','invalid type combination')
        End Select
      Class Default
        Call dash_error(0,'complex_matrix_multiply','invalid type combination')
      End Select
      M%daggered = .False.
      N%daggered = .False.
      O%daggered = .False.
    End Subroutine complex_matrix_multiply

    Subroutine real_matrix_diag(M,evals)
    Use dash_utils
      Class(real_ks_matrix), Intent(InOut) :: M
      Type(ks_vector_1D)   , Intent(InOut) :: evals
      Call bk_real_sym_diag(M%data,evals%data)
    End Subroutine real_matrix_diag

    Subroutine complex_matrix_diag(M,evals)
    Use dash_utils
      Class(complex_ks_matrix), Intent(InOut) :: M
      Type(ks_vector_1D)      , Intent(InOut) :: evals
      Call bk_comp_sym_diag(M%data,evals%data)
    End Subroutine complex_matrix_diag

    Subroutine real_matrix_copy(M,N)
    Use dash_utils
      Class(real_ks_matrix), Intent(InOut) :: M
      Class(ks_matrix)     , Intent(InOut) :: N
      Select Type(N)
      Class Is(real_ks_matrix)
        Call bk_real_matrix_copy(M%data,N%data)
      Class Default
        Call dash_error(0,'real_matrix_copy','types not coherent')
      End Select
    End Subroutine real_matrix_copy

    Subroutine complex_matrix_copy(M,N)
    Use dash_utils
      Class(complex_ks_matrix), Intent(InOut) :: M
      Class(ks_matrix)        , Intent(InOut) :: N
      Select Type(N)
      Class Is(complex_ks_matrix)
        Call bk_comp_matrix_copy(M%data,N%data)
      Class Default
        Call dash_error(0,'complex_matrix_copy','types not coherent')
      End Select
    End Subroutine complex_matrix_copy
  
    Subroutine real_matrix_compute_density(M,nocc,pg,weights,coords)
      Class(real_ks_matrix),               Intent(InOut) :: M
      Integer              ,               Intent(In   ) :: nocc
      Real(fp64)         , Dimension(:), Intent(InOut) :: pg
      Real(fp64)         , Dimension(:), Intent(In   ) :: weights
      Integer              , Dimension(3), Intent(In   ) :: coords
      Call bk_real_matrix_compute_density( M%data,nocc,pg,weights,coords )
    End Subroutine real_matrix_compute_density

    Subroutine complex_matrix_compute_density(M,nocc,pg,weights,coords)
      Class(complex_ks_matrix),               Intent(InOut) :: M
      Integer                 ,               Intent(In   ) :: nocc
      Real(fp64)            , Dimension(:), Intent(InOut) :: pg
      Real(fp64)            , Dimension(:), Intent(In   ) :: weights
      Integer                 , Dimension(3), Intent(In   ) :: coords
      Call bk_comp_matrix_compute_density( M%data,nocc,pg,weights,coords )
    End Subroutine complex_matrix_compute_density

    Subroutine real_matrix_cholesky(M)
      Class(real_ks_matrix), Intent(InOut) :: M
      Call bk_real_cholesky(M%data)
    End Subroutine real_matrix_cholesky

    Subroutine real_matrix_clean_lower(M)
      Class(real_ks_matrix), Intent(InOut) :: M
      Call bk_real_set_lower_to(M%data,0.0_fp64)
    End Subroutine real_matrix_clean_lower

    Subroutine real_matrix_invert(M)
      Class(real_ks_matrix), Intent(InOut) :: M
      Call bk_real_invert(M%data)
    End Subroutine real_matrix_invert

    Subroutine complex_matrix_cholesky(M)
      Class(complex_ks_matrix), Intent(InOut) :: M
      Call bk_comp_cholesky(M%data)
    End Subroutine complex_matrix_cholesky

    Subroutine complex_matrix_clean_lower(M)
      Class(complex_ks_matrix), Intent(InOut) :: M
      Call bk_comp_set_lower_to(M%data,Cmplx(0.0_fp64, Kind=fp64))
    End Subroutine complex_matrix_clean_lower

    Subroutine complex_matrix_invert(M)
      Class(complex_ks_matrix), Intent(InOut) :: M
      Call bk_comp_invert(M%data)
    End Subroutine complex_matrix_invert






















  !! --- HOST MEMORY STORAGE ---
    !! --- Store matrix in host memory ---
    Subroutine real_matrix_host_create(M,n)
      Class(real_ks_matrix ), Intent(InOut) :: M
      Integer               , Intent(In   ) :: n
      ! Local Variables
      Integer :: ldm
      ldm = max(1,n)
      Call bk_host_real_alloc(M%data,ldm,ldm)
    End Subroutine real_matrix_host_create

    Subroutine complex_matrix_host_create(M,n)
      Class(complex_ks_matrix), Intent(InOut) :: M
      Integer                 , Intent(In   ) :: n
      !Local Variables
      Integer :: ldm
      ldm = max(1,n)
      Call bk_host_comp_alloc(M%data,ldm,ldm)
    End Subroutine complex_matrix_host_create

    !! --- Free matrix host memory ---
    Subroutine real_matrix_host_free(M)
    !! Create the data for a real n*n matrix
      Class(real_ks_matrix ), Intent(InOut) :: M
      Call bk_host_real_dealloc(M%data)
    End Subroutine real_matrix_host_free

    Subroutine complex_matrix_host_free(M)
    !! Create the data for a complex n*n matrix
      Class(complex_ks_matrix), Intent(InOut) :: M
      Call bk_host_comp_dealloc(M%data)
    End Subroutine complex_matrix_host_free

    Subroutine real_matrix_host_copy(M,N)
    Use dash_utils
      Class(real_ks_matrix), Intent(InOut) :: M
      Class(ks_matrix)     , Intent(InOut) :: N
      Select Type(N)
      Class Is(real_ks_matrix)
        Call bk_host_real_matrix_copy(M%data,N%data)
      Class Default
        Call dash_error(0,'real_matrix_host_copy','types not coherent')
      End Select
    End Subroutine real_matrix_host_copy

    Subroutine complex_matrix_host_copy(M,N)
    Use dash_utils
      Class(complex_ks_matrix), Intent(InOut) :: M
      Class(ks_matrix)        , Intent(InOut) :: N
      Select Type(N)
      Class Is(complex_ks_matrix)
        Call bk_host_comp_matrix_copy(M%data,N%data)
      Class Default
        Call dash_error(0,'complex_matrix_host_copy','types not coherent')
      End Select
    End Subroutine complex_matrix_host_copy

    Subroutine real_matrix_host_to_dev(M)
      Class(real_ks_matrix), Intent(InOut) :: M
      Call bk_real_host_to_dev(M%data)
    End Subroutine real_matrix_host_to_dev

    Subroutine real_matrix_dev_to_host(M)
      Class(real_ks_matrix), Intent(InOut) :: M
      Call bk_real_dev_to_host(M%data)
    End Subroutine real_matrix_dev_to_host

    Subroutine complex_matrix_host_to_dev(M)
      Class(complex_ks_matrix), Intent(InOut) :: M
      Call bk_comp_host_to_dev(M%data)
    End Subroutine complex_matrix_host_to_dev

    Subroutine complex_matrix_dev_to_host(M)
      Class(complex_ks_matrix), Intent(InOut) :: M
      Call bk_comp_dev_to_host(M%data)
    End Subroutine complex_matrix_dev_to_host


    !! --- Set raw host data ---
    Subroutine real_matrix_host_set_raw_real(M,raw_data)
      Class(real_ks_matrix),                 Intent(InOut) :: M
      Real(fp64)         , Dimension(:,:), Intent(In   ) :: raw_data
      Call bk_host_set_raw_real(M%data,raw_data)
    End Subroutine real_matrix_host_set_raw_real

    Subroutine real_matrix_host_set_raw_complex(M,raw_data)
    Use dash_utils
      Class(real_ks_matrix),               Intent(InOut) :: M
      Complex(fp64)      , Dimension(:,:), Intent(In   ) :: raw_data
      Call dash_error(0,'real_matrix_host_set_raw_complex','invalid')
    End Subroutine real_matrix_host_set_raw_complex

    Subroutine complex_matrix_host_set_raw_real(M,raw_data)
    Use dash_utils
      Class(complex_ks_matrix)                , Intent(InOut) :: M
      Real(fp64)            , Dimension(:,:), Intent(In   ) :: raw_data
      Call dash_error(0,'complex_matrix_host_set_raw_real','invalid assign')
    End Subroutine complex_matrix_host_set_raw_real

    Subroutine complex_matrix_host_set_raw_complex(M,raw_data)
      Class(complex_ks_matrix)                , Intent(InOut) :: M
      Complex(fp64)         , Dimension(:,:), Intent(In   ) :: raw_data
      Call bk_host_set_raw_comp(M%data,raw_data)
    End Subroutine complex_matrix_host_set_raw_complex

    !! --- Get the raw host data ---
    Subroutine real_matrix_host_get_raw_real(M,raw_data)
      Class(real_ks_matrix),                              Intent(In   ) :: M
      Real(fp64)         , Dimension(:,:), Allocatable, Intent(InOut) :: raw_data
      Call bk_host_get_raw_real(M%data,raw_data)
    End Subroutine real_matrix_host_get_raw_real

    Subroutine real_matrix_host_get_raw_complex(M,raw_data)
    Use dash_utils
      Class(real_ks_matrix),                                Intent(In   ) :: M
      Complex(fp64)        , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call dash_error(0,'real_matrix_host_get_raw_complex','invalid')
    End Subroutine real_matrix_host_get_raw_complex

    Subroutine complex_matrix_host_get_raw_real(M,raw_data)
    Use dash_utils
      Class(complex_ks_matrix)                ,              Intent(In   ) :: M
      Real(fp64)            , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call dash_error(0,'complex_matrix_host_get_raw_complex','invalid')
    End Subroutine complex_matrix_host_get_raw_real

    Subroutine complex_matrix_host_get_raw_complex(M,raw_data)
      Class (complex_ks_matrix)                             , Intent(In   ) :: M
      Complex(fp64)          , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call bk_host_get_raw_comp(M%data,raw_data)
    End Subroutine complex_matrix_host_get_raw_complex

End Module module_matrix_v2
