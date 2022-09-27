


! #+CALL: generate_f_interface(table=qmckl_dgemm_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_dgemm &
      (context, TransA, TransB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: TransA
    character           , intent(in)  , value :: TransB
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: k
    real    (c_double ) , intent(in)  , value :: alpha
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(in)  , value :: beta
    real    (c_double ) , intent(out)         :: C(ldc,*)
    integer (c_int64_t) , intent(in)  , value :: ldc

  end function qmckl_dgemm
end interface




! #+CALL: generate_f_interface(table=qmckl_dgemm_safe_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm_safe")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_dgemm_safe &
      (context, TransA, TransB, m, n, k, alpha, A, size_A, lda, B, size_B, ldb, beta, C, size_C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: TransA
    character           , intent(in)  , value :: TransB
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: k
    real    (c_double ) , intent(in)  , value :: alpha
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: size_A
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: size_B
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(in)  , value :: beta
    real    (c_double ) , intent(out)         :: C(ldc,*)
    integer (c_int64_t) , intent(in)  , value :: size_C
    integer (c_int64_t) , intent(in)  , value :: ldc

  end function qmckl_dgemm_safe
end interface



! #+CALL: generate_f_interface(table=qmckl_adjugate_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_adjugate &
      (context, n, A, lda, B, ldb, det_l) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(out)         :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(inout)        :: det_l

  end function qmckl_adjugate
end interface



! #+CALL: generate_f_interface(table=qmckl_adjugate_safe_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate_safe")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_adjugate_safe &
      (context, n, A, size_A, lda, B, size_B, ldb, det_l) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: size_A
    real    (c_double ) , intent(out)         :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    integer (c_int64_t) , intent(in)  , value :: size_B 
    real    (c_double ) , intent(inout)        :: det_l

  end function qmckl_adjugate_safe
end interface
