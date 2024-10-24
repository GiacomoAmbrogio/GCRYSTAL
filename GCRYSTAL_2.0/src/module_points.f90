Module module_points

!===========================================================
!    @author: Giacomo Ambrogio       date: Feb 2024
!
!    DASH - Device Accelerated Solutions for HPC
!
!    Point class module
!
!    This module handles the points data structures
!    currently is just a wrapper that deals with the
!    abstract type of ks_matrix
!===========================================================

Use base_numbers
Use module_matrix_v2

  Implicit None
  Private

  Type, Public :: ks_point_1D
    Type(ks_vector_1D), Private :: vector
    Contains
    Procedure, Public :: create  => ks_point_1D_create
    Procedure, Public :: destroy => ks_point_1D_destroy
    Procedure, Public :: get_raw => ks_point_1D_get_raw
    Procedure, Public :: set_raw => ks_point_1D_set_raw
  End Type

  Type, Public :: ks_point
    Class(ks_matrix), Allocatable, Private :: matrix
    Contains
      !! --- LA Operations ---
    Procedure, Public :: cholesky    => ks_point_cholesky    !! Compute cholesky decomposition inplace, return U
    Procedure, Public :: clean_lower => ks_point_clean_lower !! Set lower-triangule to 0
    Procedure, Public :: invert      => ks_point_invert      !! Compute inverted matrix inplace (upper-triangule)
    Procedure, Public :: dagger      => ks_point_set_dagger  !! Set dagger logic based on op (see OP_? in base_numbers)
    Procedure, Public :: multiply    => ks_point_multiply    !! Multiply two matrices, self is result
    Procedure, Public :: diag        => ks_point_diag        !! Diagonalize lower-triangule inplace, requires workspace
    Procedure, Public :: copy        => ks_point_copy        !! Diagonalize lower-triangule inplace, requires workspace
      !! --- Iterative Operations ---
    Procedure, Public :: shift_diagonal => ks_point_shift_diagonal !! Perform level shift on occupied diagonal elements
      !! --- Utilities ---
    Procedure, Public :: create    => ks_point_create              !! Create the ks_point, alloc real/comp ks_matrix
    Procedure, Public :: store     => ks_point_allocate_space
    Procedure, Public :: free      => ks_point_deallocate_space
    Procedure, Public :: allocated => ks_point_allocated
    Generic  , Public :: set_raw   => set_raw_real, set_raw_comp
    Generic  , Public :: get_raw   => get_raw_real, get_raw_comp, get_raw_vect
      !! Private Implementations
    Procedure, Private :: set_raw_real => ks_point_set_raw_real
    Procedure, Private :: set_raw_comp => ks_point_set_raw_comp
    Procedure, Private :: get_raw_real => ks_point_get_raw_real
    Procedure, Private :: get_raw_vect => ks_point_get_raw_vector
    Procedure, Private :: get_raw_comp => ks_point_get_raw_comp
    Procedure, Public :: compute_tester     => ks_point_compute_tester
  End Type

  Contains
!!#############################################################################
!!        ks_point_1D type-bound procedures
!!#############################################################################
    Subroutine ks_point_1D_create(P1D,ldp)
    Use dash_utils
      Class(ks_point_1D), Intent(InOut) :: P1D
      Integer           , Intent(In   ) :: ldp
      If(ldp.lt.1) Call dash_error(0,"ks_point_1D_create","allocating vector with dim lower than 1")
      Call P1D%vector%create(ldp)
    End Subroutine ks_point_1D_create

    Subroutine ks_point_1D_destroy(P1D)
      Class(ks_point_1D), Intent(InOut) :: P1D
      Call P1D%vector%destroy()
    End Subroutine ks_point_1D_destroy

    Subroutine ks_point_1D_get_raw(P1D,raw_data)
      Class(ks_point_1D),                      Intent(InOut) :: P1D
      Real(double), Dimension(:), Allocatable, Intent(  Out) :: raw_data
      Call P1D%vector%get_raw(raw_data)
    End Subroutine ks_point_1D_get_raw

    Subroutine ks_point_1D_set_raw(P1D,raw_data)
      Class(ks_point_1D),         Intent(InOut) :: P1D
      Real(double), Dimension(:), Intent(In   ) :: raw_data
      Write(*,*) "ERROR############## ks_point_1D_set_raw NOT IMPLEMENTED"
      !Call bk_set_raw_1d(P1D%vector,raw_data)
    End Subroutine ks_point_1D_set_raw !! SYMMETRY not implmented

!!#############################################################################
!!        ks_point type-bound procedures
!!#############################################################################
!! --- LA Operations ---
    Subroutine ks_point_cholesky(P)
      Class(ks_point), Intent(InOut) :: P
      Call P%matrix%cholesky()
    End Subroutine ks_point_cholesky

    Subroutine ks_point_clean_lower(P)
      Class(ks_point), Intent(InOut) :: P
      Call P%matrix%clean_lower()
    End Subroutine ks_point_clean_lower

    Subroutine ks_point_invert(P)
      Class(ks_point), Intent(InOut) :: P
      Call P%matrix%invert()
    End Subroutine ks_point_invert

    Subroutine ks_point_set_dagger(P,op)
      Class(ks_point), Intent(InOut) :: P
      Integer        , Intent(In   ) :: op
      Call P%matrix%dagger(op)
    End Subroutine ks_point_set_dagger

    Subroutine ks_point_multiply(P,Q,R)
      Class(ks_point), Intent(InOut) :: P
      Type(ks_point) , Intent(InOut) :: Q
      Type(ks_point) , Intent(InOut) :: R
      Call P%matrix%multiply(Q%matrix,R%matrix)
    End Subroutine ks_point_multiply

    Subroutine ks_point_diag(P,evals)
      Class(ks_point)  , Intent(InOut) :: P
      Type(ks_point_1D), Intent(InOut) :: evals
      Call P%matrix%diag(evals%vector)
    End Subroutine ks_point_diag

    Subroutine ks_point_copy(P,Q)
      Class(ks_point), Intent(InOut) :: P
      Type(ks_point) , Intent(InOut) :: Q
      Call P%matrix%copy(Q%matrix)
    End Subroutine ks_point_copy

!! --- Iterative Operations ---
    Subroutine ks_point_shift_diagonal(P,shift,nocc)
      Class(ks_point), Intent(InOut) :: P
      Real(double)   , Intent(In   ) :: shift
      Integer        , Intent(In   ) :: nocc
      Call P%matrix%shift_diagonal(shift,nocc)
    End Subroutine ks_point_shift_diagonal

!! --- Utilities ---
    Subroutine ks_point_create(P,is_complex,nirr)
    Use dash_utils
      Class(ks_point), Intent(  Out) :: P
      Logical        , Intent(In   ) :: is_complex
      Integer        , Intent(In   ) :: nirr
      ! Local variables
      Integer :: ierr
      If(nirr.gt.1) Call dash_error(0,'ks_point_create','Symmetry not allowed')
      !If(Allocated(P%matrix)) Call dash_error(0,'ks_point_create','already alllocated')
      If(is_complex)Then
        Allocate(complex_ks_matrix :: P%matrix, stat=ierr)
      Else
        Allocate(real_ks_matrix    :: P%matrix, stat=ierr)
      Endif
      If(ierr.ne.0) Call dash_error(0,'ks_point_create','failed abstract allocation')
    End Subroutine ks_point_create

    Subroutine ks_point_allocate_space(P,i,ldn)
      Class(ks_point), Intent(InOut) :: P
      Integer        , Intent(In   ) :: i
      Integer        , Intent(In   ) :: ldn
      Call P%matrix%create(ldn)
    End Subroutine ks_point_allocate_space

            Subroutine ks_point_deallocate_space(P)   !! TO BE REWORKED
              Class(ks_point), Intent(InOut) :: P
              Call P%matrix%free()
            End Subroutine ks_point_deallocate_space

    Function ks_point_allocated(P) Result(res)
      Logical                     :: res
      Class(ks_point), Intent(In) :: P
      res = P%matrix%allocated()
    End Function ks_point_allocated

    Subroutine ks_point_set_raw_real(P, raw_data)
      Class(ks_point),                 Intent(InOut) :: P
      Real(double)   , Dimension(:,:), Intent(In   ) :: raw_data
      Call P%matrix%set_raw(raw_data)
    End Subroutine ks_point_set_raw_real

    Subroutine ks_point_set_raw_comp(P, raw_data)
      Class(ks_point),                 Intent(InOut) :: P
      Complex(double), Dimension(:,:), Intent(In   ) :: raw_data
      Call P%matrix%set_raw(raw_data)
    End Subroutine ks_point_set_raw_comp

    Subroutine ks_point_get_raw_real(P, raw_data)
      Class(ks_point),                              Intent(InOut) :: P
      Real(double)   , Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call P%matrix%get_raw(raw_data)
    End Subroutine ks_point_get_raw_real

    Subroutine ks_point_get_raw_vector(P,raw_vect,k_types)
      Class(ks_point),                            Intent(InOut) :: P
      Real(double)   , Dimension(:), Allocatable, Intent(  Out) :: raw_vect
      Integer        ,                            Intent(In   ) :: k_types
      Call P%matrix%get_raw(raw_vect,k_types)
    End Subroutine ks_point_get_raw_vector

    Subroutine ks_point_get_raw_comp(P, raw_data)
      Class(ks_point),                              Intent(InOut) :: P
      Complex(double), Dimension(:,:), Allocatable, Intent(  Out) :: raw_data
      Call P%matrix%get_raw(raw_data)
    End Subroutine ks_point_get_raw_comp

    Function ks_point_compute_tester(P,nocc) Result(tstr)
      Real(double)                :: tstr
      Class(ks_point), Intent(In) :: P
      Integer        , Intent(In) :: nocc
      tstr = P%matrix%compute_tester(nocc)
    End Function ks_point_compute_tester

End Module module_points



