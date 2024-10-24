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
    Type(bk_vector_1D), Private :: data
    Contains
    Procedure, Public :: create  => ks_vector_1D_create
    Procedure, Public :: destroy => ks_vector_1D_destroy
    Procedure, Public :: get_raw => ks_vector_1D_get_raw
    Procedure, Public :: set_raw => ks_vector_1D_set_raw
  End Type ks_vector_1D

  Type, Abstract, Public :: ks_matrix
  !! An abstract for the base type.
  Logical, Private :: daggered = .False.
    Contains
    ! Public methods that are overridden
    Procedure(create), Deferred, Public :: create
    Procedure(free)  , Deferred, Public :: free
    Generic  , Public :: set_raw => set_raw_real, set_raw_complex
    Generic  , Public :: get_raw => get_raw_real, get_raw_complex, get_raw_vect
    Procedure(allocated)      , Deferred, Public  :: allocated
    Procedure(diag)           , Deferred, Public  :: diag
    Procedure(copy)           , Deferred, Public  :: copy
    Procedure(compute_tester) , Deferred, Public  :: compute_tester
    Procedure(shift_diagonal) , Deferred, Public  :: shift_diagonal
    ! Private implementations
    Procedure                 ,           Public  :: dagger => ks_matrix_set_dagger ! new
    Procedure(binary_op)      , Deferred, Public  :: multiply
    Procedure(unary_op)       , Deferred, Public  :: cholesky    ! new
    Procedure(unary_op)       , Deferred, Public  :: clean_lower ! new
    Procedure(unary_op)       , Deferred, Public  :: invert      ! new
    Procedure(set_raw_real)   , Deferred, Private :: set_raw_real
    Procedure(set_raw_complex), Deferred, Private :: set_raw_complex
    Procedure(get_raw_real)   , Deferred, Private :: get_raw_real
    Procedure(get_raw_vect)   , Deferred, Private :: get_raw_vect
    Procedure(get_raw_complex), Deferred, Private :: get_raw_complex
  End type ks_matrix

  Type, Extends(ks_matrix), Public :: real_ks_matrix
  !! An instance of a distributed matrix that holds real data
    Type(bk_real_matrix), Private :: data
    Contains
    ! Public methods
    Procedure, Public :: allocated => real_matrix_allocated
    Procedure, Public :: create => real_matrix_create
    Procedure, Public :: free   => real_matrix_free
    Procedure, Public :: diag   => real_matrix_diag
    Procedure, Public :: copy   => real_matrix_copy
    Procedure, Public :: compute_tester => real_matrix_compute_tester
    Procedure, Public :: shift_diagonal => real_matrix_shift_diagonal
    ! Private implementations
    Procedure, Public  :: multiply        => real_matrix_multiply
    Procedure, Public  :: cholesky        => real_matrix_cholesky    ! new
    Procedure, Public  :: clean_lower     => real_matrix_clean_lower ! new
    Procedure, Public  :: invert          => real_matrix_invert      ! new
    Procedure, Private :: set_raw_real    => real_matrix_set_raw_real
    Procedure, Private :: set_raw_complex => real_matrix_set_raw_complex
    Procedure, Private :: get_raw_real    => real_matrix_get_raw_real
    Procedure, Private :: get_raw_vect    => real_matrix_get_raw_vector
    Procedure, Private :: get_raw_complex => real_matrix_get_raw_complex
  End type real_ks_matrix
!
  Type, Extends(ks_matrix), Public :: complex_ks_matrix
  !! An instance of a distributed matrix that holds complex data
    Type(bk_comp_matrix), Private :: data
    Contains
    ! Public methods
    Procedure, Public :: allocated => complex_matrix_allocated
    Procedure, Public :: create => complex_matrix_create                 !! Create storage for the data of the matrix 
    Procedure, Public :: free   => complex_matrix_free                 !! Create storage for the data of the matrix 
    Procedure, Public :: diag   => complex_matrix_diag
    Procedure, Public :: copy   => complex_matrix_copy
    Procedure, Public :: compute_tester => complex_matrix_compute_tester
    Procedure, Public :: shift_diagonal => complex_matrix_shift_diagonal
    ! Private implementations
    Procedure, Public :: multiply        => complex_matrix_multiply
    Procedure, Public  :: cholesky        => complex_matrix_cholesky    ! new
    Procedure, Public  :: clean_lower     => complex_matrix_clean_lower ! new
    Procedure, Public  :: invert          => complex_matrix_invert      ! new
    Procedure, Private :: set_raw_real    => complex_matrix_set_raw_real
    Procedure, Private :: set_raw_complex => complex_matrix_set_raw_complex
    Procedure, Private :: get_raw_real    => complex_matrix_get_raw_real
    Procedure, Private :: get_raw_vect    => complex_matrix_get_raw_vector
    Procedure, Private :: get_raw_complex => complex_matrix_get_raw_complex
  End type complex_ks_matrix

  Abstract Interface
    Function allocated(M) Result(res)
    Import :: ks_matrix
    Implicit None
    Logical :: res
    Class(ks_matrix), Intent(In) :: M
    End Function allocated
    Subroutine create(M,n)
    !! Create storage for the data of the matrix
    Import :: ks_matrix
    Implicit None
    Class(ks_matrix), Intent(  Out) :: M
    Integer         , Intent(In   ) :: n
    End Subroutine create
    Subroutine free(M)
    !! Create storage for the data of the matrix
    Import :: ks_matrix
    Implicit None
    Class(ks_matrix), Intent(InOut) :: M
    End Subroutine free
    Function compute_tester(M,nocc) Result(tstr)
    Import :: double
    Import :: ks_matrix
    Real(double) :: tstr
    Class(ks_matrix), Intent(In) :: M
    Integer         , Intent(In) :: nocc
    End Function compute_tester
    Subroutine shift_diagonal(M,shift,nocc)
    Import :: double
    Import :: ks_matrix
    Class(ks_matrix), Intent(InOut) :: M
    Real(double)    , Intent(In) :: shift
    Integer         , Intent(In) :: nocc
    End Subroutine shift_diagonal
    Subroutine get_raw_real(M,raw_data)
    !! Get raw data from a real array
    Import :: double
    Import :: ks_matrix
    Implicit None
    Class(ks_matrix),                              Intent(In   ) :: M
    Real(double)    , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
    End Subroutine get_raw_real
    Subroutine get_raw_vect(M,raw_vect,kt)
    !! Get raw data from a real array
    Import :: double
    Import :: ks_matrix
    Implicit None
    Class(ks_matrix),                            Intent(In   ) :: M
    Real(double)    , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
    Integer         ,                            Intent(In   ) :: kt
    End Subroutine get_raw_vect
    Subroutine get_raw_complex(M,raw_data)
    !! Get raw data from a complex array
    Import :: double
    Import :: ks_matrix
    Implicit None
    Class(ks_matrix),                              Intent(In   ) :: M
    Complex(double) , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
    End Subroutine get_raw_complex
    Subroutine set_raw_real(M,raw_data)
    !! Set raw data from a real array
    Import :: double
    Import :: ks_matrix
    Implicit None
    Class(ks_matrix),                 Intent(InOut) :: M
    Real(double)    , Dimension(:,:), Intent(In   ) :: raw_data
    End Subroutine set_raw_real
    Subroutine set_raw_complex(M,raw_data)
    !! Set raw data from a compex array
    Import :: double
    Import :: ks_matrix
    Implicit None
    Class(ks_matrix),                 Intent(InOut) :: M
    Complex(double) , Dimension(:,:), Intent(In   ) :: raw_data
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
!!#############################################################################
!!        ks_vector_1D type-bound procedures
!!#############################################################################
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
      Real(double), Dimension(:), Allocatable, Intent(  Out) :: raw_data
      Call bk_get_raw_1d(V1D%data,raw_data)
    End Subroutine ks_vector_1D_get_raw

    Subroutine ks_vector_1D_set_raw(V1D,raw_data)
      Class(ks_vector_1D),         Intent(InOut) :: V1D
      Real(double), Dimension(:), Intent(In   ) :: raw_data
      Write(*,*) "ERROR############## ks_vector_1D_set_raw NOT IMPLEMENTED"
      !Call bk_set_raw_1d(P1D%vector,raw_data)
    End Subroutine ks_vector_1D_set_raw

    Subroutine ks_matrix_set_dagger(M,op)
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
  Logical :: res
  Class(real_ks_matrix), Intent(In) :: M
    res = bk_real_allocated(M%data)
  End Function real_matrix_allocated

  Function complex_matrix_allocated(M) Result(res)
  Logical :: res
  Class(complex_ks_matrix), Intent(In) :: M
    res = bk_comp_allocated(M%data)
  End Function complex_matrix_allocated

  !! Over-ridding routines
    !! --- Store matrix in memory ---
    Subroutine real_matrix_create(M,n)
    !! Create the data for a real n*n matrix
      Class(real_ks_matrix ), Intent(  Out) :: M
      Integer               , Intent(In   ) :: n
      ! Local Variables
      Integer :: ldm
      ldm = max(1,n)
      Call bk_real_alloc(M%data,ldm,ldm)
    End Subroutine real_matrix_create
    Subroutine complex_matrix_create(M,n)
    !! Create the data for a complex n*n matrix
      Class(complex_ks_matrix), Intent(  Out) :: M
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
      Real(double)         , Dimension(:,:), Intent(In   ) :: raw_data
      Call bk_set_raw_real(M%data,raw_data)
    End Subroutine real_matrix_set_raw_real
    Subroutine real_matrix_set_raw_complex(M,raw_data)
    !! Set the raw data for M
    Use dash_utils
      Class(real_ks_matrix),               Intent(InOut) :: M
      Complex(double)      , Dimension(:,:), Intent(In   ) :: raw_data
      Call dash_error(0,'real_matrix_set_raw_complex','invalid')
    End Subroutine real_matrix_set_raw_complex
    Subroutine complex_matrix_set_raw_real(M,raw_data)
    !! set the raw data for M
    Use dash_utils
      Class(complex_ks_matrix)                , Intent(InOut) :: M
      Real(double)            , Dimension(:,:), Intent(In   ) :: raw_data
      Call dash_error(0,'complex_matrix_set_raw_real','invalid assign')
      ! It is possible to assign a real matrix to a complex one,
      ! but still it is not the cleanest way to do it
      !M%data = raw_data
    End Subroutine complex_matrix_set_raw_real
    Subroutine complex_matrix_set_raw_complex(M,raw_data)
    !! Set the raw data for M
      Class(complex_ks_matrix)                , Intent(InOut) :: M
      Complex(double)         , Dimension(:,:), Intent(In   ) :: raw_data
      Call bk_set_raw_comp(M%data,raw_data)
    End Subroutine complex_matrix_set_raw_complex

    !! --- Get the raw data ---
    Subroutine real_matrix_get_raw_real(M,raw_data)
    !! Get the raw data for M
      Class(real_ks_matrix),                              Intent(In   ) :: M
      Real(double)         , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call bk_get_raw_real(M%data,raw_data)
    End Subroutine real_matrix_get_raw_real
    Subroutine real_matrix_get_raw_vector(M,raw_vect,kt)
    !! get raw data from M, put in CRYSTAL like 1D vector, std data ordering
      Class(real_ks_matrix),                            Intent(In   ) :: M
      Real(double)         , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
      Integer              ,                            Intent(In   ) :: kt
      Call bk_get_raw_real_vect(M%data,raw_vect,kt)
    End Subroutine real_matrix_get_raw_vector
    Subroutine real_matrix_get_raw_complex(M,raw_data)
    !! Get the raw data for A
    Use dash_utils
      Class(real_ks_matrix),                                Intent(In   ) :: M
      Complex(double)        , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call dash_error(0,'real_matrix_get_raw_complex','invalid')
    End Subroutine real_matrix_get_raw_complex
    Subroutine complex_matrix_get_raw_real(M,raw_data)
    !! Get the raw data for M
    Use dash_utils
      Class(complex_ks_matrix)                ,              Intent(In   ) :: M
      Real(double)            , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call dash_error(0,'complex_matrix_get_raw_complex','invalid')
    End Subroutine complex_matrix_get_raw_real
    Subroutine complex_matrix_get_raw_vector(M,raw_vect,kt)
    !! Get the raw data for M
    Use dash_utils
      Class(complex_ks_matrix),                            Intent(In   ) :: M
      Real(double)            , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
      Integer                 ,                            Intent(In   ) :: kt
      ! Local Variables
      Integer, Parameter :: ONLY_REAL_KS = 1
      Integer, Parameter :: BOTH_TYPES_KS = 2
      Integer            :: i, j
      Call bk_get_raw_comp_vect(M%data,raw_vect,kt)
    End Subroutine complex_matrix_get_raw_vector
    Subroutine complex_matrix_get_raw_complex(M,raw_data)
    !! Get the raw data for M
      Class (complex_ks_matrix)                             , Intent(In   ) :: M
      Complex(double)          , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call bk_get_raw_comp(M%data,raw_data)
    End Subroutine complex_matrix_get_raw_complex

    Function real_matrix_compute_tester(M,nocc) Result(tstr)
      Real(double)                      :: tstr
      Class(real_ks_matrix), Intent(In) :: M
      Integer              , Intent(In) :: nocc

      Call bk_compute_real_tester_sub(M%data,nocc,tstr)
    End Function real_matrix_compute_tester

    Function complex_matrix_compute_tester(M,nocc) Result(tstr)
      Real(double)                         :: tstr
      Class(complex_ks_matrix), Intent(In) :: M
      Integer                 , Intent(In) :: nocc

      Call bk_compute_complex_tester_sub(M%data,nocc,tstr)
    End Function complex_matrix_compute_tester

    Subroutine real_matrix_shift_diagonal(M,shift,nocc)
      Class(real_ks_matrix), Intent(InOut) :: M
      Real(double)         , Intent(In   ) :: shift
      Integer              , Intent(In   ) :: nocc
      Call bk_real_shift_diag(M%data,shift,nocc)
    End Subroutine real_matrix_shift_diagonal

    Subroutine complex_matrix_shift_diagonal(M,shift,nocc)
      Class(complex_ks_matrix), Intent(InOut) :: M
      Real(double)            , Intent(In   ) :: shift
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

    !! --- Choleski Routines ---
    Subroutine real_matrix_cholesky(M)
      Class(real_ks_matrix), Intent(InOut) :: M
      Call bk_real_cholesky(M%data)
    End Subroutine real_matrix_cholesky

    Subroutine real_matrix_clean_lower(M)
      Class(real_ks_matrix), Intent(InOut) :: M
      Call bk_real_set_lower_to(M%data,0.0_double)
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
      Call bk_comp_set_lower_to(M%data,Cmplx(0.0_double, Kind=double))
    End Subroutine complex_matrix_clean_lower

    Subroutine complex_matrix_invert(M)
      Class(complex_ks_matrix), Intent(InOut) :: M
      Call bk_comp_invert(M%data)
    End Subroutine complex_matrix_invert

End Module module_matrix_v2
