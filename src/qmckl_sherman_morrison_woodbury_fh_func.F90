! Fortran interface                                               :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_sherman_morrison
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_sherman_morrison_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sherman_morrison &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)

    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    integer (c_int64_t) , intent(in) , value :: N_updates
    real    (c_double ) , intent(in)         :: Updates(N_updates*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(N_updates)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_sherman_morrison
end interface

! Fortran interface                                              :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_woodbury_2
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_woodbury_2_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_2 &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    real    (c_double ) , intent(in)         :: Updates(2*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(2)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_woodbury_2
end interface

! Fortran interface                                               :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_woodbury_3
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_woodbury_3_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_3 &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    real    (c_double ) , intent(in)         :: Updates(3*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(3)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_woodbury_3
end interface

! Fortran interface                                              :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_sherman_morrison_splitting
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_sherman_morrison_splitting_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sherman_morrison_splitting &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    integer (c_int64_t) , intent(in) , value :: N_updates
    real    (c_double ) , intent(in)         :: Updates(N_updates*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(N_updates)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_sherman_morrison_splitting
end interface

! Fortran interface                                               :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_sherman_morrison_smw32s
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_sherman_morrison_smw32s_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sherman_morrison_smw32s &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    integer (c_int64_t) , intent(in) , value :: N_updates
    real    (c_double ) , intent(in)         :: Updates(N_updates*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(N_updates)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_sherman_morrison_smw32s
end interface
