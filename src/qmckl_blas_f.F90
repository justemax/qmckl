integer function qmckl_dgemm_f(context, TransA, TransB, &
     m, n, k, alpha, A, LDA, B, LDB, beta, C, LDC) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  character             , intent(in)  :: TransA, TransB
  integer*8             , intent(in)  :: m, n, k
  double precision      , intent(in)  :: alpha, beta
  integer*8             , intent(in)  :: lda
  double precision      , intent(in)  :: A(lda,*)
  integer*8             , intent(in)  :: ldb
  double precision      , intent(in)  :: B(ldb,*)
  integer*8             , intent(in)  :: ldc
  double precision      , intent(out) :: C(ldc,*)

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (m <= 0_8) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (n <= 0_8) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  if (k <= 0_8) then
     info = QMCKL_INVALID_ARG_6
     return
  endif

  if (LDA <= 0) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  if (LDB <= 0) then
     info = QMCKL_INVALID_ARG_11
     return
  endif

  if (LDC <= 0) then
     info = QMCKL_INVALID_ARG_14
     return
  endif

  call dgemm(transA, transB, int(m,4), int(n,4), int(k,4), &
       alpha, A, int(LDA,4), B, int(LDB,4), beta, C, int(LDC,4))

end function qmckl_dgemm_f

! C interface                                                    :noexport:

!     #+CALL: generate_c_interface(table=qmckl_dgemm_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm")

!     #+RESULTS:

integer(c_int32_t) function qmckl_dgemm &
    (context, TransA, TransB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
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

  integer(c_int32_t), external :: qmckl_dgemm_f
  info = qmckl_dgemm_f &
         (context, TransA, TransB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)

end function qmckl_dgemm

integer function qmckl_dgemm_safe_f(context, TransA, TransB, &
     m, n, k, alpha, A, size_A, LDA, B, size_B, LDB, beta, C, size_C, LDC) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  character             , intent(in)  :: TransA, TransB
  integer*8             , intent(in)  :: m, n, k
  double precision      , intent(in)  :: alpha, beta
  integer*8             , intent(in)  :: lda
  integer*8             , intent(in)  :: size_A
  double precision      , intent(in)  :: A(lda,*)
  integer*8             , intent(in)  :: ldb
  integer*8             , intent(in)  :: size_B
  double precision      , intent(in)  :: B(ldb,*)
  integer*8             , intent(in)  :: ldc
  integer*8             , intent(in)  :: size_C
  double precision      , intent(out) :: C(ldc,*)

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (m <= 0_8) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (n <= 0_8) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  if (k <= 0_8) then
     info = QMCKL_INVALID_ARG_6
     return
  endif

  if (LDA <= 0) then
     info = QMCKL_INVALID_ARG_10
     return
  endif

  if (LDB <= 0) then
     info = QMCKL_INVALID_ARG_13
     return
  endif

  if (LDC <= 0) then
     info = QMCKL_INVALID_ARG_17
     return
  endif

  if (size_A <= 0) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  if (size_B <= 0) then
     info = QMCKL_INVALID_ARG_12
     return
  endif

  if (size_C <= 0) then
     info = QMCKL_INVALID_ARG_16
     return
  endif

  call dgemm(transA, transB, int(m,4), int(n,4), int(k,4), &
       alpha, A, int(LDA,4), B, int(LDB,4), beta, C, int(LDC,4))

end function qmckl_dgemm_safe_f

! C interface                                                    :noexport:

!     #+CALL: generate_c_interface(table=qmckl_dgemm_safe_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm_safe")

!     #+RESULTS:

integer(c_int32_t) function qmckl_dgemm_safe &
    (context, TransA, TransB, m, n, k, alpha, A, size_A, lda, B, size_B, ldb, beta, C, size_C, ldc) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
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

  integer(c_int32_t), external :: qmckl_dgemm_safe_f
  info = qmckl_dgemm_safe_f &
         (context, TransA, TransB, m, n, k, alpha, A, size_A, lda, B, size_B, ldb, beta, C, size_C, ldc)

end function qmckl_dgemm_safe

integer function qmckl_adjugate_f(context, na, A, LDA, B, ldb, det_l) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  double precision, intent(in)          :: A (LDA,*)
  integer*8, intent(in)                 :: LDA
  double precision, intent(out)         :: B (LDB,*)
  integer*8, intent(in)                 :: LDB
  integer*8, intent(in)                 :: na
  double precision, intent(inout)       :: det_l

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (na <= 0_8) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (LDA <= 0_8) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (LDA < na) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  select case (na)
  case (5)
     call adjugate5(A,LDA,B,LDB,na,det_l)
  case (4)
     call adjugate4(A,LDA,B,LDB,na,det_l)
  case (3)
     call adjugate3(A,LDA,B,LDB,na,det_l)
  case (2)
     call adjugate2(A,LDA,B,LDB,na,det_l)
  case (1)
    det_l = a(1,1)
    b(1,1) = 1.d0
  case default
     call adjugate_general(context, na, A, LDA, B, LDB, det_l)
  end select

end function qmckl_adjugate_f

subroutine adjugate2(A,LDA,B,LDB,na,det_l)
  implicit none
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDA,na)
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(inout) :: det_l

  double precision :: C(2,2)

  call cofactor2(A,LDA,C,2_8,na,det_l)

  B(1,1) = C(1,1)
  B(2,1) = C(1,2)
  B(1,2) = C(2,1)
  B(2,2) = C(2,2)

end subroutine adjugate2

subroutine adjugate3(a,LDA,B,LDB,na,det_l)
  implicit none
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDA,na)
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(inout) :: det_l

  double precision :: C(4,3)

  call cofactor3(A,LDA,C,4_8,na,det_l)

  B(1,1) = C(1,1)
  B(1,2) = C(2,1)
  B(1,3) = C(3,1)
  B(2,1) = C(1,2)
  B(2,2) = C(2,2)
  B(2,3) = C(3,2)
  B(3,1) = C(1,3)
  B(3,2) = C(2,3)
  B(3,3) = C(3,3)

end subroutine adjugate3

subroutine adjugate4(a,LDA,B,LDB,na,det_l)
  implicit none
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDA,na)
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(inout) :: det_l

  double precision :: C(4,4)

  call cofactor4(A,LDA,C,4_8,na,det_l)
  B(1,1) = C(1,1)
  B(1,2) = C(2,1)
  B(1,3) = C(3,1)
  B(1,4) = C(4,1)
  B(2,1) = C(1,2)
  B(2,2) = C(2,2)
  B(2,3) = C(3,2)
  B(2,4) = C(4,2)
  B(3,1) = C(1,3)
  B(3,2) = C(2,3)
  B(3,3) = C(3,3)
  B(3,4) = C(4,3)
  B(4,1) = C(1,4)
  B(4,2) = C(2,4)
  B(4,3) = C(3,4)
  B(4,4) = C(4,4)

end subroutine adjugate4

subroutine adjugate5(A,LDA,B,LDB,na,det_l)
  implicit none
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDA,na)
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(inout) :: det_l

  double precision  :: C(8,5)

  call cofactor5(A,LDA,C,8_8,na,det_l)

  B(1,1) = C(1,1)
  B(1,2) = C(2,1)
  B(1,3) = C(3,1)
  B(1,4) = C(4,1)
  B(1,5) = C(5,1)
  B(2,1) = C(1,2)
  B(2,2) = C(2,2)
  B(2,3) = C(3,2)
  B(2,4) = C(4,2)
  B(2,5) = C(5,2)
  B(3,1) = C(1,3)
  B(3,2) = C(2,3)
  B(3,3) = C(3,3)
  B(3,4) = C(4,3)
  B(3,5) = C(5,3)
  B(4,1) = C(1,4)
  B(4,2) = C(2,4)
  B(4,3) = C(3,4)
  B(4,4) = C(4,4)
  B(4,5) = C(5,4)
  B(5,1) = C(1,5)
  B(5,2) = C(2,5)
  B(5,3) = C(3,5)
  B(5,4) = C(4,5)
  B(5,5) = C(5,5)

end subroutine adjugate5

subroutine cofactor2(a,LDA,b,LDB,na,det_l)
  implicit none
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDA,na)
  integer*8, intent(in)           :: LDA, LDB
  integer*8        :: na
  double precision :: det_l

  det_l  =  a(1,1)*a(2,2) - a(1,2)*a(2,1)
  b(1,1) =  a(2,2)
  b(2,1) = -a(2,1)
  b(1,2) = -a(1,2)
  b(2,2) =  a(1,1)
end subroutine cofactor2

subroutine cofactor3(a,LDA,b,LDB,na,det_l)
  implicit none
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDA,na)
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(inout) :: det_l
  integer :: i

  det_l = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
         -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
         +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

  b(1,1) =  a(2,2)*a(3,3) - a(2,3)*a(3,2)
  b(2,1) =  a(2,3)*a(3,1) - a(2,1)*a(3,3)
  b(3,1) =  a(2,1)*a(3,2) - a(2,2)*a(3,1)

  b(1,2) =  a(1,3)*a(3,2) - a(1,2)*a(3,3)
  b(2,2) =  a(1,1)*a(3,3) - a(1,3)*a(3,1)
  b(3,2) =  a(1,2)*a(3,1) - a(1,1)*a(3,2)

  b(1,3) =  a(1,2)*a(2,3) - a(1,3)*a(2,2)
  b(2,3) =  a(1,3)*a(2,1) - a(1,1)*a(2,3)
  b(3,3) =  a(1,1)*a(2,2) - a(1,2)*a(2,1)

end subroutine cofactor3

subroutine cofactor4(a,LDA,b,LDB,na,det_l)
  implicit none
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDA,na)
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(inout) :: det_l
  integer :: i,j
  det_l =  a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
                  -a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
                  +a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) &
          -a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
                  -a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
                  +a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))) &
          +a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
                  -a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
                  +a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))) &
          -a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))  &
                  -a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))  &
                  +a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

  b(1,1) =  a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))
  b(2,1) = -a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))
  b(3,1) =  a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))
  b(4,1) = -a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))-a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))

  b(1,2) = -a(1,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(1,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(1,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))
  b(2,2) =  a(1,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(1,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(1,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))
  b(3,2) = -a(1,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(1,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(1,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))
  b(4,2) =  a(1,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(1,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(1,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))

  b(1,3) =  a(1,2)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))-a(1,3)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))+a(1,4)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))
  b(2,3) = -a(1,1)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))+a(1,3)*(a(2,1)*a(4,4)-a(2,4)*a(4,1))-a(1,4)*(a(2,1)*a(4,3)-a(2,3)*a(4,1))
  b(3,3) =  a(1,1)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))-a(1,2)*(a(2,1)*a(4,4)-a(2,4)*a(4,1))+a(1,4)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))
  b(4,3) = -a(1,1)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))+a(1,2)*(a(2,1)*a(4,3)-a(2,3)*a(4,1))-a(1,3)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))

  b(1,4) = -a(1,2)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))+a(1,3)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))-a(1,4)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
  b(2,4) =  a(1,1)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))-a(1,3)*(a(2,1)*a(3,4)-a(2,4)*a(3,1))+a(1,4)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
  b(3,4) = -a(1,1)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))+a(1,2)*(a(2,1)*a(3,4)-a(2,4)*a(3,1))-a(1,4)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  b(4,4) =  a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

end subroutine cofactor4

subroutine cofactor5(A,LDA,B,LDB,na,det_l)
  implicit none
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDA,na)
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(inout) :: det_l
  integer :: i,j

 det_l = a(1,1)*(a(2,2)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*( &
 a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))- &
 a(2,3)*(a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)- &
 a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(2,4)*(a(3,2)*( &
 a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+ &
 a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,5)*(a(3,2)*(a(4,3)*a(5,4)- &
 a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)* &
 a(5,3)-a(4,3)*a(5,2))))-a(1,2)*(a(2,1)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)* &
 a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)- &
 a(4,4)*a(5,3)))-a(2,3)*(a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*( &
 a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+ &
 a(2,4)*(a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)- &
 a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(2,5)*(a(3,1)*( &
 a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+ &
 a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))+a(1,3)*(a(2,1)*(a(3,2)*(a(4,4)* &
 a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*( &
 a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(2,2)*(a(3,1)*(a(4,4)*a(5,5)-a(4,5)* &
 a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)- &
 a(4,4)*a(5,1)))+a(2,4)*(a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*( &
 a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))- &
 a(2,5)*(a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)- &
 a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))-a(1,4)*(a(2,1)*( &
 a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)* &
 a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)*(a(3,1)*(a(4,3)* &
 a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*( &
 a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)*(a(3,1)*(a(4,2)*a(5,5)-a(4,5)* &
 a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)- &
 a(4,2)*a(5,1)))-a(2,5)*(a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*( &
 a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))+ &
 a(1,5)*(a(2,1)*(a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)* &
 a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)*( &
 a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)* &
 a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)*(a(3,1)*(a(4,2)* &
 a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*( &
 a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(2,4)*(a(3,1)*(a(4,2)*a(5,3)-a(4,3)* &
 a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)- &
 a(4,2)*a(5,1))))

 b(1,1) = &
 (a(2,2)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))-a(2,3)* &
 (a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(2,4)* &
 (a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,5)* &
 (a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))))
 b(2,1) = &
 (-a(2,1)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))+a(2,3)* &
 (a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))-a(2,4)* &
 (a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,5)* &
 (a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))
 b(3,1) = &
 (a(2,1)*(a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(2,2)* &
 (a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+a(2,4)* &
 (a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(2,5)* &
 (a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(4,1) = &
 (-a(2,1)*(a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))+a(2,2)* &
 (a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(2,3)* &
 (a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))+a(2,5)* &
 (a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(5,1) = &
 (a(2,1)*(a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)* &
 (a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)* &
 (a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(2,4)* &
 (a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))

 b(1,2) = &
 (-a(1,2)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))+a(1,3)* &
 (a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(1,4)* &
 (a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))+a(1,5)* &
 (a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))))
 b(2,2) = &
 (a(1,1)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))-a(1,3)* &
 (a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+a(1,4)* &
 (a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(1,5)* &
 (a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))
 b(3,2) = &
 (-a(1,1)*(a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(1,2)* &
 (a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))-a(1,4)* &
 (a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))+a(1,5)* &
 (a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(4,2) = &
 (a(1,1)*(a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(1,2)* &
 (a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(1,3)* &
 (a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(1,5)* &
 (a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(5,2) = &
 (-a(1,1)*(a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))+a(1,2)* &
 (a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(1,3)* &
 (a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))+a(1,4)* &
 (a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))

 b(1,3) = &
 (a(1,2)*(a(2,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(2,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))-a(1,3)* &
 (a(2,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(2,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(1,4)* &
 (a(2,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(2,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(2,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(1,5)* &
 (a(2,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(2,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(2,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))))
 b(2,3) = &
 (-a(1,1)*(a(2,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(2,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))+a(1,3)* &
 (a(2,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))-a(1,4)* &
 (a(2,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(2,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(1,5)* &
 (a(2,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(2,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(2,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))
 b(3,3) = &
 (a(1,1)*(a(2,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(2,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(1,2)* &
 (a(2,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+a(1,4)* &
 (a(2,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(2,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(1,5)* &
 (a(2,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(2,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(2,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(4,3) = &
 (-a(1,1)*(a(2,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(2,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(2,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))+a(1,2)* &
 (a(2,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(2,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(1,3)* &
 (a(2,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(2,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))+a(1,5)* &
 (a(2,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(2,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(2,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(5,3) = &
 (a(1,1)*(a(2,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(2,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(2,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(1,2)* &
 (a(2,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(2,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(2,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(1,3)* &
 (a(2,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(2,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(2,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(1,4)* &
 (a(2,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(2,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(2,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))

 b(1,4) = &
 (-a(1,2)*(a(2,3)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))+a(2,5)*(a(3,3)*a(5,4)-a(3,4)*a(5,3)))+a(1,3)* &
 (a(2,2)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))+a(2,5)*(a(3,2)*a(5,4)-a(3,4)*a(5,2)))-a(1,4)* &
 (a(2,2)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))-a(2,3)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))+a(2,5)*(a(3,2)*a(5,3)-a(3,3)*a(5,2)))+a(1,5)* &
 (a(2,2)*(a(3,3)*a(5,4)-a(3,4)*a(5,3))-a(2,3)*(a(3,2)*a(5,4)-a(3,4)*a(5,2))+a(2,4)*(a(3,2)*a(5,3)-a(3,3)*a(5,2))))
 b(2,4) = &
 (a(1,1)*(a(2,3)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))+a(2,5)*(a(3,3)*a(5,4)-a(3,4)*a(5,3)))-a(1,3)* &
 (a(2,1)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,4)-a(3,4)*a(5,1)))+a(1,4)* &
 (a(2,1)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))-a(2,3)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,3)-a(3,3)*a(5,1)))-a(1,5)* &
 (a(2,1)*(a(3,3)*a(5,4)-a(3,4)*a(5,3))-a(2,3)*(a(3,1)*a(5,4)-a(3,4)*a(5,1))+a(2,4)*(a(3,1)*a(5,3)-a(3,3)*a(5,1))))
 b(3,4) = &
 (-a(1,1)*(a(2,2)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))+a(2,5)*(a(3,2)*a(5,4)-a(3,4)*a(5,2)))+a(1,2)* &
 (a(2,1)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,4)-a(3,4)*a(5,1)))-a(1,4)* &
 (a(2,1)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))-a(2,2)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,2)-a(3,2)*a(5,1)))+a(1,5)* &
 (a(2,1)*(a(3,2)*a(5,4)-a(3,4)*a(5,2))-a(2,2)*(a(3,1)*a(5,4)-a(3,4)*a(5,1))+a(2,4)*(a(3,1)*a(5,2)-a(3,2)*a(5,1))))
 b(4,4) = &
 (a(1,1)*(a(2,2)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))-a(2,3)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))+a(2,5)*(a(3,2)*a(5,3)-a(3,3)*a(5,2)))-a(1,2)* &
 (a(2,1)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))-a(2,3)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,3)-a(3,3)*a(5,1)))+a(1,3)* &
 (a(2,1)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))-a(2,2)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,2)-a(3,2)*a(5,1)))-a(1,5)* &
 (a(2,1)*(a(3,2)*a(5,3)-a(3,3)*a(5,2))-a(2,2)*(a(3,1)*a(5,3)-a(3,3)*a(5,1))+a(2,3)*(a(3,1)*a(5,2)-a(3,2)*a(5,1))))
 b(5,4) = &
 (-a(1,1)*(a(2,2)*(a(3,3)*a(5,4)-a(3,4)*a(5,3))-a(2,3)*(a(3,2)*a(5,4)-a(3,4)*a(5,2))+a(2,4)*(a(3,2)*a(5,3)-a(3,3)*a(5,2)))+a(1,2)* &
 (a(2,1)*(a(3,3)*a(5,4)-a(3,4)*a(5,3))-a(2,3)*(a(3,1)*a(5,4)-a(3,4)*a(5,1))+a(2,4)*(a(3,1)*a(5,3)-a(3,3)*a(5,1)))-a(1,3)* &
 (a(2,1)*(a(3,2)*a(5,4)-a(3,4)*a(5,2))-a(2,2)*(a(3,1)*a(5,4)-a(3,4)*a(5,1))+a(2,4)*(a(3,1)*a(5,2)-a(3,2)*a(5,1)))+a(1,4)* &
 (a(2,1)*(a(3,2)*a(5,3)-a(3,3)*a(5,2))-a(2,2)*(a(3,1)*a(5,3)-a(3,3)*a(5,1))+a(2,3)*(a(3,1)*a(5,2)-a(3,2)*a(5,1))))

 b(1,5) = &
 (a(1,2)*(a(2,3)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))+a(2,5)*(a(3,3)*a(4,4)-a(3,4)*a(4,3)))-a(1,3)* &
 (a(2,2)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))+a(2,5)*(a(3,2)*a(4,4)-a(3,4)*a(4,2)))+a(1,4)* &
 (a(2,2)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))-a(2,3)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))+a(2,5)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))-a(1,5)* &
 (a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))))
 b(2,5) = &
 (-a(1,1)*(a(2,3)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))+a(2,5)*(a(3,3)*a(4,4)-a(3,4)*a(4,3)))+a(1,3)* &
 (a(2,1)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,4)-a(3,4)*a(4,1)))-a(1,4)* &
 (a(2,1)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))-a(2,3)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))+a(1,5)* &
 (a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))))
 b(3,5) = &
 (a(1,1)*(a(2,2)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))+a(2,5)*(a(3,2)*a(4,4)-a(3,4)*a(4,2)))-a(1,2)* &
 (a(2,1)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,4)-a(3,4)*a(4,1)))+a(1,4)* &
 (a(2,1)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))-a(2,2)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))-a(1,5)* &
 (a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))))
 b(4,5) = &
 (-a(1,1)*(a(2,2)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))-a(2,3)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))+a(2,5)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))+a(1,2)* &
 (a(2,1)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))-a(2,3)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))-a(1,3)* &
 (a(2,1)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))-a(2,2)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))+a(1,5)* &
 (a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))))
 b(5,5) = &
 (a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))-a(1,2)* &
 (a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))+a(1,3)* &
 (a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))-a(1,4)* &
 (a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))))

end



! #+CALL: generate_c_interface(table=qmckl_adjugate_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate")

! #+RESULTS:

integer(c_int32_t) function qmckl_adjugate &
    (context, n, A, lda, B, ldb, det_l) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: n
  real    (c_double ) , intent(in)          :: A(lda,*)
  integer (c_int64_t) , intent(in)  , value :: lda
  real    (c_double ) , intent(out)         :: B(ldb,*)
  integer (c_int64_t) , intent(in)  , value :: ldb
  real    (c_double ) , intent(inout)        :: det_l

  integer(c_int32_t), external :: qmckl_adjugate_f
  info = qmckl_adjugate_f &
         (context, n, A, lda, B, ldb, det_l)

end function qmckl_adjugate

subroutine adjugate_general(context, na, A, LDA, B, LDB, det_l)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in) :: context
  double precision, intent(in)       :: A (LDA,na)
  integer*8, intent(in)              :: LDA
  double precision, intent(out)      :: B (LDB,na)
  integer*8, intent(in)              :: LDB
  integer*8, intent(in)              :: na
  double precision, intent(inout)    :: det_l

  double precision :: work(LDA*max(na,64))
  integer          :: inf
  integer          :: ipiv(LDA)
  integer          :: lwork
  integer(8)       :: i, j

B(1:na,1:na) = A(1:na,1:na)

call dgetrf(na, na, B, LDB, ipiv, inf )

det_l = 1.d0
j=0_8
do i=1,na
 j = j+min(abs(ipiv(i)-i),1)
 det_l = det_l*B(i,i)
enddo

if (iand(j,1_8) /= 0_8)  then
  det_l = -det_l
endif

lwork = SIZE(work)
call dgetri(na, B, LDB, ipiv, work, lwork, inf )

B(:,:) = B(:,:)*det_l

end subroutine adjugate_general

integer function qmckl_adjugate_safe_f(context, &
     na, A, size_A, LDA, B, size_B, LDB, det_l) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  double precision, intent(in)          :: A (LDA,*)
  integer*8, intent(in)                 :: size_A 
  integer*8, intent(in)                 :: LDA
  double precision, intent(out)         :: B (LDB,*)
  integer*8, intent(in)                 :: size_B 
  integer*8, intent(in)                 :: LDB
  integer*8, intent(in)                 :: na
  double precision, intent(inout)       :: det_l

  integer, external :: qmckl_adjugate_f

  info = QMCKL_SUCCESS

  if (size_A < na) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (size_B <= 0_8) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  info = qmckl_adjugate_f(context, na, A, LDA, B, LDB, det_l) 

  if (info == QMCKL_INVALID_ARG_4) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  if (info == QMCKL_INVALID_ARG_6) then
     info = QMCKL_INVALID_ARG_8
     return
  endif

end function qmckl_adjugate_safe_f

! C interface

!     #+CALL: generate_c_interface(table=qmckl_adjugate_safe_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate_safe")

!     #+RESULTS:

integer(c_int32_t) function qmckl_adjugate_safe &
    (context, n, A, size_A, lda, B, size_B, ldb, det_l) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
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

  integer(c_int32_t), external :: qmckl_adjugate_safe_f
  info = qmckl_adjugate_safe_f &
         (context, n, A, size_A, lda, B, size_B, ldb, det_l)

end function qmckl_adjugate_safe
