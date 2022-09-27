

! #+CALL: generate_f_interface(table=qmckl_distance_sq_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_distance_sq &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transa
    character           , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(out)         :: C(ldc,n)
    integer (c_int64_t) , intent(in)  , value :: ldc

  end function qmckl_distance_sq
end interface



! #+CALL: generate_f_interface(table=qmckl_distance_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_distance &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transa
    character           , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(out)         :: C(ldc,n)
    integer (c_int64_t) , intent(in)  , value :: ldc

  end function qmckl_distance
end interface



! #+CALL: generate_f_interface(table=qmckl_distance_rescaled_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_distance_rescaled &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc, rescale_factor_kappa) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transa
    character           , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(out)         :: C(ldc,n)
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)  , value :: rescale_factor_kappa

  end function qmckl_distance_rescaled
end interface



! #+CALL: generate_f_interface(table=qmckl_distance_rescaled_deriv_e_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_distance_rescaled_deriv_e &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc, rescale_factor_kappa) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transa
    character           , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(out)         :: C(ldc,n,4)
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)  , value :: rescale_factor_kappa

  end function qmckl_distance_rescaled_deriv_e
end interface
