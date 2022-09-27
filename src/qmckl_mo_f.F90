! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_mo_basis_mo_value
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_mo_basis_mo_value_args
!     | Variable        | Type                        | In/Out | Description                                     |
!     |-----------------+-----------------------------+--------+-------------------------------------------------|
!     | ~context~       | ~qmckl_context~             | in     | Global state                                    |
!     | ~ao_num~        | ~int64_t~                   | in     | Number of AOs                                   |
!     | ~mo_num~        | ~int64_t~                   | in     | Number of MOs                                   |
!     | ~point_num~     | ~int64_t~                   | in     | Number of points                                |
!     | ~coefficient_t~ | ~double[mo_num][ao_num]~    | in     | Transpose of the AO to MO transformation matrix |
!     | ~ao_value~      | ~double[point_num][ao_num]~ | in     | Value of the AOs                                |
!     | ~mo_value~      | ~double[point_num][mo_num]~ | out    | Value of the MOs                                |


!     The matrix of AO values is very sparse, so we use a sparse-dense
!     matrix multiplication instead of a dgemm, as exposed in
!     https://dx.doi.org/10.1007/978-3-642-38718-0_14.




integer function qmckl_compute_mo_basis_mo_value_doc_f(context, &
     ao_num, mo_num, point_num, &
     coefficient_t, ao_value, mo_value) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: ao_num, mo_num
  integer*8             , intent(in)  :: point_num
  double precision      , intent(in)  :: ao_value(ao_num,point_num)
  double precision      , intent(in)  :: coefficient_t(mo_num,ao_num)
  double precision      , intent(out) :: mo_value(mo_num,point_num)
  integer*8 :: i,j,k
  double precision :: c1, c2, c3, c4, c5

  integer*8 :: LDA, LDB, LDC

  info = QMCKL_SUCCESS
  if (.True.)  then    ! fast algorithm
     do j=1,point_num
        mo_value(:,j) = 0.d0
        do k=1,ao_num
           c1 = ao_value(k,j)
           if (c1 /= 0.d0) then
              do i=1,mo_num
                 mo_value(i,j) = mo_value(i,j) + coefficient_t(i,k) * c1
              end do
           end if
        end do
     end do
     
  else ! dgemm for checking

    LDA = size(coefficient_t,1)
    LDB = size(ao_value,1) 
    LDC = size(mo_value,1)

    info = qmckl_dgemm(context,'N', 'N', mo_num, point_num, ao_num, 1.d0,     &
                                    coefficient_t, LDA, ao_value, LDB, &
                                    0.d0, mo_value, LDC)

  end if

end function qmckl_compute_mo_basis_mo_value_doc_f



! #+CALL: generate_c_interface(table=qmckl_mo_basis_mo_value_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value_doc"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_mo_basis_mo_value_doc &
    (context, ao_num, mo_num, point_num, coefficient_t, ao_value, mo_value) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: mo_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  real    (c_double ) , intent(in)          :: coefficient_t(ao_num,mo_num)
  real    (c_double ) , intent(in)          :: ao_value(ao_num,point_num)
  real    (c_double ) , intent(out)         :: mo_value(mo_num,point_num)

  integer(c_int32_t), external :: qmckl_compute_mo_basis_mo_value_doc_f
  info = qmckl_compute_mo_basis_mo_value_doc_f &
         (context, ao_num, mo_num, point_num, coefficient_t, ao_value, mo_value)

end function qmckl_compute_mo_basis_mo_value_doc

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_mo_basis_mo_vgl
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_mo_basis_mo_vgl_args
!     | Variable            | Type                           | In/Out | Description                                     |
!     |---------------------+--------------------------------+--------+-------------------------------------------------|
!     | ~context~           | ~qmckl_context~                | in     | Global state                                    |
!     | ~ao_num~            | ~int64_t~                      | in     | Number of AOs                                   |
!     | ~mo_num~            | ~int64_t~                      | in     | Number of MOs                                   |
!     | ~point_num~         | ~int64_t~                      | in     | Number of points                                |
!     | ~coefficient_t~     | ~double[mo_num][ao_num]~       | in     | Transpose of the AO to MO transformation matrix |
!     | ~ao_vgl~            | ~double[point_num][5][ao_num]~ | in     | Value, gradients and Laplacian of the AOs       |
!     | ~mo_vgl~            | ~double[point_num][5][mo_num]~ | out    | Value, gradients and Laplacian of the MOs       |


!     The matrix of AO values is very sparse, so we use a sparse-dense
!     matrix multiplication instead of a dgemm, as exposed in
!     https://dx.doi.org/10.1007/978-3-642-38718-0_14.




integer function qmckl_compute_mo_basis_mo_vgl_doc_f(context, &
     ao_num, mo_num, point_num, &
     coefficient_t, ao_vgl, mo_vgl) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: ao_num, mo_num
  integer*8             , intent(in)  :: point_num
  double precision      , intent(in)  :: ao_vgl(ao_num,5,point_num)
  double precision      , intent(in)  :: coefficient_t(mo_num,ao_num)
  double precision      , intent(out) :: mo_vgl(mo_num,5,point_num)
  integer*8 :: i,j,k
  double precision :: c1, c2, c3, c4, c5

  do j=1,point_num
     mo_vgl(:,:,j) = 0.d0
     do k=1,ao_num
        if (ao_vgl(k,1,j) /= 0.d0) then
           c1 = ao_vgl(k,1,j)
           c2 = ao_vgl(k,2,j)
           c3 = ao_vgl(k,3,j)
           c4 = ao_vgl(k,4,j)
           c5 = ao_vgl(k,5,j)
           do i=1,mo_num
              mo_vgl(i,1,j) = mo_vgl(i,1,j) + coefficient_t(i,k) * c1
              mo_vgl(i,2,j) = mo_vgl(i,2,j) + coefficient_t(i,k) * c2
              mo_vgl(i,3,j) = mo_vgl(i,3,j) + coefficient_t(i,k) * c3
              mo_vgl(i,4,j) = mo_vgl(i,4,j) + coefficient_t(i,k) * c4
              mo_vgl(i,5,j) = mo_vgl(i,5,j) + coefficient_t(i,k) * c5
           end do
        end if
     end do
  end do
  info = QMCKL_SUCCESS

! info = qmckl_dgemm(context,'N', 'N', mo_num, point_num, ao_num, 1.d0, &
!      coefficient_t, int(size(coefficient_t,1),8),      &
!      ao_vgl, int(size(ao_vgl,1),8), 0.d0,                  &
!      mo_vgl, int(size(mo_vgl,1),8))

end function qmckl_compute_mo_basis_mo_vgl_doc_f



! #+CALL: generate_c_interface(table=qmckl_mo_basis_mo_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl_doc"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_mo_basis_mo_vgl_doc &
    (context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: mo_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  real    (c_double ) , intent(in)          :: coefficient_t(ao_num,mo_num)
  real    (c_double ) , intent(in)          :: ao_vgl(ao_num,5,point_num)
  real    (c_double ) , intent(out)         :: mo_vgl(mo_num,5,point_num)

  integer(c_int32_t), external :: qmckl_compute_mo_basis_mo_vgl_doc_f
  info = qmckl_compute_mo_basis_mo_vgl_doc_f &
         (context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl)

end function qmckl_compute_mo_basis_mo_vgl_doc
