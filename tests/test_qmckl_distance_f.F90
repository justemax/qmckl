integer(qmckl_exit_code) function test_qmckl_distance_sq(context) bind(C)

  use qmckl
  use qmckl_verificarlo_f
  use iso_c_binding

  implicit none

  integer(qmckl_context), intent(in), value :: context
  logical(C_BOOL) :: vfc_err

  double precision, allocatable :: A(:,:), B(:,:), C(:,:)
  integer*8                     :: m, n, LDA, LDB, LDC
  double precision              :: x
  integer*8                     :: i,j

  m = 5
  n = 6
  LDA = m
  LDB = n
  LDC = 5

  allocate( A(LDA,m), B(LDB,n), C(LDC,n) )
  do j=1,m
     do i=1,m
        A(i,j) = -10.d0 + dble(i+j)
     end do
  end do
  do j=1,n
     do i=1,n
        B(i,j) = -1.d0 + dble(i*j)
     end do
  end do

  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'X', 't', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe("distance"//C_NULL_CHAR, "distance_sq_Xt_2_2"//C_NULL_CHAR, DBLE(C(2,2)))

  if (test_qmckl_distance_sq == 0) return


  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 't', 'X', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe("distance"//C_NULL_CHAR, "distance_sq_tX_2_2"//C_NULL_CHAR, DBLE(C(2,2)))

  if (test_qmckl_distance_sq == 0) return


  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'T', 't', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance"//C_NULL_CHAR, "distance_sq_Tt_2_2"//C_NULL_CHAR, DBLE(C(2,2)), DBLE(0), DBLE(1.d-14))

  if (test_qmckl_distance_sq == 0) return


  test_qmckl_distance_sq = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  (A(i,1)-B(j,1))**2 + &
             (A(i,2)-B(j,2))**2 + &
             (A(i,3)-B(j,3))**2
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'n', 'T', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance"//C_NULL_CHAR, "distance_sq_nT_2_2"//C_NULL_CHAR, DBLE(C(2,2)), DBLE(0), DBLE(1.d-14))


  test_qmckl_distance_sq = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  (A(1,i)-B(j,1))**2 + &
             (A(2,i)-B(j,2))**2 + &
             (A(3,i)-B(j,3))**2
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'T', 'n', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err =  qmckl_probe_check("distance"//C_NULL_CHAR, "distance_sq_Tn_2_2"//C_NULL_CHAR, DBLE(C(2,2)), DBLE(0), DBLE(1.d-14))

  if (test_qmckl_distance_sq == 0) return

  test_qmckl_distance_sq = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  (A(i,1)-B(1,j))**2 + &
             (A(i,2)-B(2,j))**2 + &
             (A(i,3)-B(3,j))**2
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'n', 'N', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance"//C_NULL_CHAR, "distance_sq_nN_2_2"//C_NULL_CHAR, DBLE(C(2,2)), DBLE(0), DBLE(1.d-14))

  test_qmckl_distance_sq = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  (A(1,i)-B(1,j))**2 + &
             (A(2,i)-B(2,j))**2 + &
             (A(3,i)-B(3,j))**2
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14 ) return
     end do
  end do

  test_qmckl_distance_sq = QMCKL_SUCCESS

  deallocate(A,B,C)
end function test_qmckl_distance_sq

integer(qmckl_exit_code) function test_qmckl_dist(context) bind(C)

  use qmckl
  use qmckl_verificarlo_f
  use iso_c_binding

  implicit none

  integer(qmckl_context), intent(in), value :: context
  logical(C_BOOL) :: vfc_err

  double precision, allocatable :: A(:,:), B(:,:), C(:,:)
  integer*8                     :: m, n, LDA, LDB, LDC
  double precision              :: x
  integer*8                     :: i,j

  m = 5
  n = 6
  LDA = m
  LDB = n
  LDC = 5

  allocate( A(LDA,m), B(LDB,n), C(LDC,n) )

  do j=1,m
     do i=1,m
        A(i,j) = -10.d0 + dble(i+j)
     end do
  end do
  do j=1,n
     do i=1,n
        B(i,j) = -1.d0 + dble(i*j)
     end do
  end do

  test_qmckl_dist = &
       qmckl_distance(context, 'X', 't', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe("distance"//C_NULL_CHAR, "distance_Xt_2_2"//C_NULL_CHAR, DBLE(C(2,2)))

  if (test_qmckl_dist == 0) return

  test_qmckl_dist = &
       qmckl_distance(context, 't', 'X', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe("distance"//C_NULL_CHAR, "distance_tX_2_2"//C_NULL_CHAR, DBLE(C(2,2)))

  if (test_qmckl_dist == 0) return

  test_qmckl_dist = &
       qmckl_distance(context, 'T', 't', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance"//C_NULL_CHAR, "distance_Tt_2_2"//C_NULL_CHAR, DBLE(C(2,2)), DBLE(0), DBLE(1.d-14))

  if (test_qmckl_dist == 0) return


  test_qmckl_dist = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  dsqrt((A(i,1)-B(j,1))**2 + &
                   (A(i,2)-B(j,2))**2 + &
                   (A(i,3)-B(j,3))**2)
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_dist = &
       qmckl_distance(context, 'n', 'T', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance"//C_NULL_CHAR, "distance_nT_2_2"//C_NULL_CHAR, DBLE(C(2,2)), DBLE(0), DBLE(1.d-14))

  if (test_qmckl_dist == 0) return


  test_qmckl_dist = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x = dsqrt((A(1,i)-B(j,1))**2 + &
                  (A(2,i)-B(j,2))**2 + &
                  (A(3,i)-B(j,3))**2)
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_dist = &
       qmckl_distance(context, 'T', 'n', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance"//C_NULL_CHAR, "distance_Tn_2_2"//C_NULL_CHAR, DBLE(C(2,2)), DBLE(0), DBLE(1.d-14))

  if (test_qmckl_dist == 0) return


  test_qmckl_dist = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  dsqrt((A(i,1)-B(1,j))**2 + &
                   (A(i,2)-B(2,j))**2 + &
                   (A(i,3)-B(3,j))**2)
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_dist = &
       qmckl_distance(context, 'n', 'N', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance"//C_NULL_CHAR, "distance_nN_2_2"//C_NULL_CHAR, DBLE(C(2,2)), DBLE(0), DBLE(1.d-14))

  if (test_qmckl_dist == 0) return


  test_qmckl_dist = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x = dsqrt((A(1,i)-B(1,j))**2 + &
                  (A(2,i)-B(2,j))**2 + &
                  (A(3,i)-B(3,j))**2)
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_dist = QMCKL_SUCCESS

  deallocate(A,B,C)
end function test_qmckl_dist
