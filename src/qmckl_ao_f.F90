#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

integer function qmckl_ao_gaussian_vgl_f(context, X, R, n, A, VGL, ldv) result(info)
  use qmckl
  implicit none
  integer*8 , intent(in)  :: context
  real*8    , intent(in)  :: X(3), R(3)
  integer*8 , intent(in)  :: n
  real*8    , intent(in)  :: A(n)
  real*8    , intent(out) :: VGL(ldv,5)
  integer*8 , intent(in)  :: ldv

  integer*8         :: i,j
  real*8            :: Y(3), r2, t, u, v

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (n <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (ldv < n) then
     info = QMCKL_INVALID_ARG_7
     return
  endif


  do i=1,3
     Y(i) = X(i) - R(i)
  end do
  r2 = Y(1)*Y(1) + Y(2)*Y(2) + Y(3)*Y(3)

  do i=1,n
     VGL(i,1) = dexp(-A(i) * r2)
  end do

  do i=1,n
     VGL(i,5) = A(i) * VGL(i,1)
  end do

  t = -2.d0 * ( X(1) - R(1) )
  u = -2.d0 * ( X(2) - R(2) )
  v = -2.d0 * ( X(3) - R(3) )

  do i=1,n
     VGL(i,2) = t * VGL(i,5)
     VGL(i,3) = u * VGL(i,5)
     VGL(i,4) = v * VGL(i,5)
  end do

  t = 4.d0 * r2
  do i=1,n
     VGL(i,5) = (t * A(i) - 6.d0) *  VGL(i,5)
  end do

end function qmckl_ao_gaussian_vgl_f

integer(c_int32_t) function qmckl_ao_gaussian_vgl(context, X, R, n, A, VGL, ldv) &
     bind(C) result(info)
  use, intrinsic :: iso_c_binding
  implicit none
  integer (c_int64_t) , intent(in) , value :: context
  real    (c_double)  , intent(in)         :: X(3), R(3)
  integer (c_int64_t) , intent(in) , value :: n
  real    (c_double)  , intent(in)         :: A(n)
  real    (c_double)  , intent(out)        :: VGL(ldv,5)
  integer (c_int64_t) , intent(in) , value :: ldv

  integer, external :: qmckl_ao_gaussian_vgl_f
  info = qmckl_ao_gaussian_vgl_f(context, X, R, n, A, VGL, ldv)
end function qmckl_ao_gaussian_vgl

integer function qmckl_compute_ao_basis_primitive_gaussian_vgl_f( &
     context, prim_num, point_num, nucl_num,                       &
     nucleus_prim_index, coord, nucl_coord,                  &
     expo, primitive_vgl)                                         &
     result(info)

  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: prim_num
  integer*8             , intent(in)  :: nucl_num
  integer*8             , intent(in)  :: point_num
  integer*8             , intent(in)  :: nucleus_prim_index(nucl_num+1)
  double precision      , intent(in)  :: coord(point_num,3)
  double precision      , intent(in)  :: nucl_coord(nucl_num,3)
  double precision      , intent(in)  :: expo(prim_num)
  double precision      , intent(out) :: primitive_vgl(prim_num,5,point_num)

  integer*8 :: inucl, iprim, ipoint
  double precision :: x, y, z, two_a, ar2, r2, v, cutoff

  info = QMCKL_SUCCESS

  ! Don't compute exponentials when the result will be almost zero.
  cutoff = 27.631021115928547  ! -dlog(1.d-12)

  do inucl=1,nucl_num
     ! C is zero-based, so shift bounds by one
     do iprim = nucleus_prim_index(inucl)+1, nucleus_prim_index(inucl+1)
        do ipoint = 1, point_num
           x = coord(ipoint,1) - nucl_coord(inucl,1)
           y = coord(ipoint,2) - nucl_coord(inucl,2)
           z = coord(ipoint,3) - nucl_coord(inucl,3)

           r2 = x*x + y*y + z*z
           ar2 = expo(iprim)*r2
           if (ar2 > cutoff) cycle

           v = dexp(-ar2)
           two_a = -2.d0 * expo(iprim) * v

           primitive_vgl(iprim, 1, ipoint) = v
           primitive_vgl(iprim, 2, ipoint) = two_a * x
           primitive_vgl(iprim, 3, ipoint) = two_a * y
           primitive_vgl(iprim, 4, ipoint) = two_a * z
           primitive_vgl(iprim, 5, ipoint) = two_a * (3.d0 - 2.d0*ar2)

        end do
     end do
  end do

end function qmckl_compute_ao_basis_primitive_gaussian_vgl_f



! #+CALL: generate_c_interface(table=qmckl_ao_basis_primitive_gaussian_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_basis_primitive_gaussian_vgl"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_ao_basis_primitive_gaussian_vgl &
    (context, &
     prim_num, &
     point_num, &
     nucl_num, &
     nucleus_prim_index, &
     coord, &
     nucl_coord, &
     expo, &
     primitive_vgl) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: prim_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)          :: nucleus_prim_index(nucl_num)
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(in)          :: expo(prim_num)
  real    (c_double ) , intent(out)         :: primitive_vgl(prim_num,5,point_num)

  integer(c_int32_t), external :: qmckl_compute_ao_basis_primitive_gaussian_vgl_f
  info = qmckl_compute_ao_basis_primitive_gaussian_vgl_f &
         (context, &
     prim_num, &
     point_num, &
     nucl_num, &
     nucleus_prim_index, &
     coord, &
     nucl_coord, &
     expo, &
     primitive_vgl)

end function qmckl_compute_ao_basis_primitive_gaussian_vgl

integer function qmckl_compute_ao_basis_shell_gaussian_vgl_f( &
     context, prim_num, shell_num, point_num, nucl_num,       &
     nucleus_shell_num, nucleus_index, nucleus_range,         &
     shell_prim_index, shell_prim_num, coord, nucl_coord,     &
     expo, coef_normalized, shell_vgl)                        &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: prim_num
  integer*8             , intent(in)  :: shell_num
  integer*8             , intent(in)  :: nucl_num
  integer*8             , intent(in)  :: point_num
  integer*8             , intent(in)  :: nucleus_shell_num(nucl_num)
  integer*8             , intent(in)  :: nucleus_index(nucl_num)
  double precision      , intent(in)  :: nucleus_range(nucl_num)
  integer*8             , intent(in)  :: shell_prim_index(shell_num)
  integer*8             , intent(in)  :: shell_prim_num(shell_num)
  double precision      , intent(in)  :: coord(point_num,3)
  double precision      , intent(in)  :: nucl_coord(nucl_num,3)
  double precision      , intent(in)  :: expo(prim_num)
  double precision      , intent(in)  :: coef_normalized(prim_num)
  double precision      , intent(out) :: shell_vgl(shell_num,5,point_num)

  integer*8 :: inucl, iprim, ipoint, ishell
  integer*8 :: ishell_start, ishell_end
  integer*8 :: iprim_start , iprim_end
  double precision :: x, y, z, two_a, ar2, r2, v, cutoff

  info = QMCKL_SUCCESS

  ! Don't compute exponentials when the result will be almost zero.
  ! TODO : Use numerical precision here
  cutoff = 27.631021115928547  !-dlog(1.d-12)

  do ipoint = 1, point_num

     do inucl=1,nucl_num

        x = coord(ipoint,1) - nucl_coord(inucl,1)
        y = coord(ipoint,2) - nucl_coord(inucl,2)
        z = coord(ipoint,3) - nucl_coord(inucl,3)

        r2 = x*x + y*y + z*z

        if (r2 > cutoff*nucleus_range(inucl)) then
           cycle
        end if

        ! C is zero-based, so shift bounds by one
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)

        do ishell=ishell_start, ishell_end

           shell_vgl(ishell, 1, ipoint) = 0.d0
           shell_vgl(ishell, 2, ipoint) = 0.d0
           shell_vgl(ishell, 3, ipoint) = 0.d0
           shell_vgl(ishell, 4, ipoint) = 0.d0
           shell_vgl(ishell, 5, ipoint) = 0.d0

           iprim_start = shell_prim_index(ishell) + 1
           iprim_end   = shell_prim_index(ishell) + shell_prim_num(ishell)

           do iprim = iprim_start, iprim_end

              ar2 = expo(iprim)*r2
              if (ar2 > cutoff) then
                 cycle
              end if

              v = coef_normalized(iprim) * dexp(-ar2)
              two_a = -2.d0 * expo(iprim) * v

              shell_vgl(ishell, 1, ipoint) = &
                   shell_vgl(ishell, 1, ipoint) + v

              shell_vgl(ishell, 2, ipoint) = &
                   shell_vgl(ishell, 2, ipoint) + two_a * x

              shell_vgl(ishell, 3, ipoint) = &
                   shell_vgl(ishell, 3, ipoint) + two_a * y

              shell_vgl(ishell, 4, ipoint) = &
                   shell_vgl(ishell, 4, ipoint) + two_a * z

              shell_vgl(ishell, 5, ipoint) = &
                   shell_vgl(ishell, 5, ipoint) + two_a * (3.d0 - 2.d0*ar2)

           end do

        end do
     end do

  end do

end function qmckl_compute_ao_basis_shell_gaussian_vgl_f



! #+CALL: generate_c_interface(table=qmckl_ao_basis_shell_gaussian_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_basis_shell_gaussian_vgl"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_ao_basis_shell_gaussian_vgl &
    (context, &
     prim_num, &
     shell_num, &
     point_num, &
     nucl_num, &
     nucleus_shell_num, &
     nucleus_index, &
     nucleus_range, &
     shell_prim_index, &
     shell_prim_num, &
     coord, &
     nucl_coord, &
     expo, &
     coef_normalized, &
     shell_vgl) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: prim_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_index(shell_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_num(shell_num)
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(in)          :: expo(prim_num)
  real    (c_double ) , intent(in)          :: coef_normalized(prim_num)
  real    (c_double ) , intent(out)         :: shell_vgl(shell_num,5,point_num)

  integer(c_int32_t), external :: qmckl_compute_ao_basis_shell_gaussian_vgl_f
  info = qmckl_compute_ao_basis_shell_gaussian_vgl_f &
         (context, &
     prim_num, &
     shell_num, &
     point_num, &
     nucl_num, &
     nucleus_shell_num, &
     nucleus_index, &
     nucleus_range, &
     shell_prim_index, &
     shell_prim_num, &
     coord, &
     nucl_coord, &
     expo, &
     coef_normalized, &
     shell_vgl)

end function qmckl_compute_ao_basis_shell_gaussian_vgl

integer function qmckl_ao_power_f(context, n, X, LMAX, P, ldp) result(info)
  use qmckl
  implicit none
  integer*8 , intent(in)  :: context
  integer*8 , intent(in)  :: n
  real*8    , intent(in)  :: X(n)
  integer   , intent(in)  :: LMAX(n)
  real*8    , intent(out) :: P(ldp,n)
  integer*8 , intent(in)  :: ldp

  integer*8  :: i,k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (n <= ldp) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  k = MAXVAL(LMAX)
  if (LDP < k) then
     info = QMCKL_INVALID_ARG_6
     return
  endif

  if (k <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  do i=1,n
     P(1,i) = X(i)
     do k=2,LMAX(i)
        P(k,i) = P(k-1,i) * X(i)
     end do
  end do

end function qmckl_ao_power_f



! #+CALL: generate_c_interface(table=qmckl_ao_power_args,rettyp=get_value("CRetType"),fname="qmckl_ao_power")

! #+RESULTS:

integer(c_int32_t) function qmckl_ao_power &
    (context, n, X, LMAX, P, ldp) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: n
  real    (c_double ) , intent(in)          :: X(n)
  integer (c_int32_t) , intent(in)          :: LMAX(n)
  real    (c_double ) , intent(out)         :: P(ldp,n)
  integer (c_int64_t) , intent(in)  , value :: ldp

  integer(c_int32_t), external :: qmckl_ao_power_f
  info = qmckl_ao_power_f &
         (context, n, X, LMAX, P, ldp)

end function qmckl_ao_power

integer function qmckl_ao_polynomial_vgl_doc_f (context,    &
     X, R, lmax, n, L, ldl, VGL, ldv) result(info)
  use qmckl
  implicit none
  integer*8 , intent(in)  :: context
  real*8    , intent(in)  :: X(3), R(3)
  integer   , intent(in)  :: lmax
  integer*8 , intent(out) :: n
  integer   , intent(out) :: L(ldl,(lmax+1)*(lmax+2)*(lmax+3)/6)
  integer*8 , intent(in)  :: ldl
  real*8    , intent(out) :: VGL(ldv,(lmax+1)*(lmax+2)*(lmax+3)/6)
  integer*8 , intent(in)  :: ldv

  integer*8         :: i,j
  integer           :: a,b,c,d
  real*8            :: Y(3)
  real*8            :: pows(-2:lmax,3)
  double precision  :: xy, yz, xz
  double precision  :: da, db, dc, dd

  info = 0

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (lmax < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (ldl < 3) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (ldv < 5) then
     info = QMCKL_INVALID_ARG_9
     return
  endif


  do i=1,3
     Y(i) = X(i) - R(i)
  end do

  if (lmax == 0) then
     VGL(1,1) = 1.d0
     VGL(2:5,1) = 0.d0
     l(1:3,1) = 0
     n=1
  else if (lmax > 0) then
     pows(-2:0,1:3) = 1.d0
     do i=1,lmax
        pows(i,1) = pows(i-1,1) * Y(1)
        pows(i,2) = pows(i-1,2) * Y(2)
        pows(i,3) = pows(i-1,3) * Y(3)
     end do

     VGL(1:5,1:4) = 0.d0
     l  (1:3,1:4) = 0

     VGL(1  ,1  ) = 1.d0
     VGL(1:5,2:4) = 0.d0

     l  (1,2) = 1
     VGL(1,2) = pows(1,1)
     VGL(2,2) = 1.d0

     l  (2,3) = 1
     VGL(1,3) = pows(1,2)
     VGL(3,3) = 1.d0

     l  (3,4) = 1
     VGL(1,4) = pows(1,3)
     VGL(4,4) = 1.d0

     n=4
  endif

  ! l>=2
  dd = 2.d0
  do d=2,lmax
     da = dd
     do a=d,0,-1
        db = dd-da
        do b=d-a,0,-1
           c  = d  - a  - b
           dc = dd - da - db
           n = n+1

           l(1,n) = a
           l(2,n) = b
           l(3,n) = c

           xy = pows(a,1) * pows(b,2)
           yz = pows(b,2) * pows(c,3)
           xz = pows(a,1) * pows(c,3)

           VGL(1,n) = xy * pows(c,3)

           xy = dc * xy
           xz = db * xz
           yz = da * yz

           VGL(2,n) = pows(a-1,1) * yz
           VGL(3,n) = pows(b-1,2) * xz
           VGL(4,n) = pows(c-1,3) * xy

           VGL(5,n) = &
                (da-1.d0) * pows(a-2,1) * yz + &
                (db-1.d0) * pows(b-2,2) * xz + &
                (dc-1.d0) * pows(c-2,3) * xy

           db = db - 1.d0
        end do
        da = da - 1.d0
     end do
     dd = dd + 1.d0
  end do

  info = QMCKL_SUCCESS

end function qmckl_ao_polynomial_vgl_doc_f



! #+CALL: generate_c_interface(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_vgl_doc" )

! #+RESULTS:

integer(c_int32_t) function qmckl_ao_polynomial_vgl_doc &
    (context, X, R, lmax, n, L, ldl, VGL, ldv) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  real    (c_double ) , intent(in)          :: X(3)
  real    (c_double ) , intent(in)          :: R(3)
  integer (c_int32_t) , intent(in)  , value :: lmax
  integer (c_int64_t) , intent(inout)        :: n
  integer (c_int32_t) , intent(out)         :: L(ldl,n)
  integer (c_int64_t) , intent(in)  , value :: ldl
  real    (c_double ) , intent(out)         :: VGL(ldv,n)
  integer (c_int64_t) , intent(in)  , value :: ldv

  integer(c_int32_t), external :: qmckl_ao_polynomial_vgl_doc_f
  info = qmckl_ao_polynomial_vgl_doc_f &
         (context, X, R, lmax, n, L, ldl, VGL, ldv)

end function qmckl_ao_polynomial_vgl_doc

integer function qmckl_ao_polynomial_transp_vgl_doc_f (context,    &
     X, R, lmax, n, L, ldl, VGL, ldv) result(info)
  use qmckl
  implicit none
  integer*8 , intent(in)  :: context
  real*8    , intent(in)  :: X(3), R(3)
  integer   , intent(in)  :: lmax
  integer*8 , intent(out) :: n
  integer   , intent(out) :: L(ldl,(lmax+1)*(lmax+2)*(lmax+3)/6)
  integer*8 , intent(in)  :: ldl
  real*8    , intent(out) :: VGL(ldv,5)
  integer*8 , intent(in)  :: ldv

  integer*8         :: i,j
  integer           :: a,b,c,d
  real*8            :: Y(3)
  real*8            :: pows(-2:21,3) ! lmax < 22
  double precision  :: xy, yz, xz
  double precision  :: da, db, dc, dd

  info = 0

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (lmax < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (ldl < 3) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (ldv < (lmax+1)*(lmax+2)*(lmax+3)/6) then
     info = QMCKL_INVALID_ARG_9
     return
  endif


  if (lmax > 0) then

     do i=1,3
        Y(i) = X(i) - R(i)
     end do
     pows(-2:0,1:3) = 1.d0
     do i=1,lmax
        pows(i,1) = pows(i-1,1) * Y(1)
        pows(i,2) = pows(i-1,2) * Y(2)
        pows(i,3) = pows(i-1,3) * Y(3)
     end do

     l  (1:3,1:4) = 0
     VGL(1:4,1:5) = 0.d0

     VGL(1  ,1  ) = 1.d0

     l  (1,2) = 1
     VGL(2,1) = Y(1)
     VGL(2,2) = 1.d0

     l  (2,3) = 1
     VGL(3,1) = Y(2)
     VGL(3,3) = 1.d0

     l  (3,4) = 1
     VGL(4,1) = Y(3)
     VGL(4,4) = 1.d0

     n=4
  else
     VGL(1,1) = 1.d0
     VGL(1,2:5) = 0.d0
     l(1:3,1) = 0
     n=1
     return
  endif

  ! l>=2
  dd = 2.d0
  do d=2,lmax
     da = dd
     do a=d,0,-1
        db = dd-da
        do b=d-a,0,-1
           c  = d  - a  - b
           dc = dd - da - db
           n = n+1

           xy = pows(a,1) * pows(b,2)
           yz = pows(b,2) * pows(c,3)
           xz = pows(a,1) * pows(c,3)

           l(1,n) = a
           l(2,n) = b
           l(3,n) = c

           VGL(n,1) = xy * pows(c,3)

           xy = dc * xy
           xz = db * xz
           yz = da * yz

           VGL(n,2) = pows(a-1,1) * yz
           VGL(n,3) = pows(b-1,2) * xz
           VGL(n,4) = pows(c-1,3) * xy

           VGL(n,5) = &
                (da-1.d0) * pows(a-2,1) * yz + &
                (db-1.d0) * pows(b-2,2) * xz + &
                (dc-1.d0) * pows(c-2,3) * xy

           db = db - 1.d0
        end do
        da = da - 1.d0
     end do
     dd = dd + 1.d0
  end do

  info = QMCKL_SUCCESS

end function qmckl_ao_polynomial_transp_vgl_doc_f



! #+CALL: generate_c_interface(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_transp_vgl_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_ao_polynomial_transp_vgl_doc &
    (context, X, R, lmax, n, L, ldl, VGL, ldv) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  real    (c_double ) , intent(in)          :: X(3)
  real    (c_double ) , intent(in)          :: R(3)
  integer (c_int32_t) , intent(in)  , value :: lmax
  integer (c_int64_t) , intent(inout)        :: n
  integer (c_int32_t) , intent(out)         :: L(ldl,n)
  integer (c_int64_t) , intent(in)  , value :: ldl
  real    (c_double ) , intent(out)         :: VGL(ldv,n)
  integer (c_int64_t) , intent(in)  , value :: ldv

  integer(c_int32_t), external :: qmckl_ao_polynomial_transp_vgl_doc_f
  info = qmckl_ao_polynomial_transp_vgl_doc_f &
         (context, X, R, lmax, n, L, ldl, VGL, ldv)

end function qmckl_ao_polynomial_transp_vgl_doc

! Unoptimized version
!      #+NAME: qmckl_ao_value_args_doc
!     | Variable              | Type                              | In/Out | Description                                  |
!     |-----------------------+-----------------------------------+--------+----------------------------------------------|
!     | ~context~             | ~qmckl_context~                   | in     | Global state                                 |
!     | ~ao_num~              | ~int64_t~                         | in     | Number of AOs                                |
!     | ~shell_num~           | ~int64_t~                         | in     | Number of shells                             |
!     | ~point_num~           | ~int64_t~                         | in     | Number of points                             |
!     | ~nucl_num~            | ~int64_t~                         | in     | Number of nuclei                             |
!     | ~coord~               | ~double[3][point_num]~            | in     | Coordinates                                  |
!     | ~nucl_coord~          | ~double[3][nucl_num]~             | in     | Nuclear  coordinates                         |
!     | ~nucleus_index~       | ~int64_t[nucl_num]~               | in     | Index of the 1st shell of each nucleus       |
!     | ~nucleus_shell_num~   | ~int64_t[nucl_num]~               | in     | Number of shells per nucleus                 |
!     | ~nucleus_range~       | ~double[nucl_num]~                | in     | Range beyond which all is zero               |
!     | ~nucleus_max_ang_mom~ | ~int32_t[nucl_num]~               | in     | Maximum angular momentum per nucleus         |
!     | ~shell_ang_mom~       | ~int32_t[shell_num]~              | in     | Angular momentum of each shell               |
!     | ~ao_factor~           | ~double[ao_num]~                  | in     | Normalization factor of the AOs              |
!     | ~shell_vgl~           | ~double[point_num][5][shell_num]~ | in     | Value, gradients and Laplacian of the shells |
!     | ~ao_value~            | ~double[point_num][ao_num]~       | out    | Values of the AOs                            |


integer function qmckl_compute_ao_value_doc_f(context, &
     ao_num, shell_num, point_num, nucl_num, &
     coord, nucl_coord, nucleus_index, nucleus_shell_num, &
     nucleus_range, nucleus_max_ang_mom, shell_ang_mom, &
     ao_factor, shell_vgl, ao_value) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: ao_num
  integer*8             , intent(in)  :: shell_num
  integer*8             , intent(in)  :: point_num
  integer*8             , intent(in)  :: nucl_num
  double precision      , intent(in)  :: coord(point_num,3)
  double precision      , intent(in)  :: nucl_coord(nucl_num,3)
  integer*8             , intent(in)  :: nucleus_index(nucl_num)
  integer*8             , intent(in)  :: nucleus_shell_num(nucl_num)
  double precision      , intent(in)  :: nucleus_range(nucl_num)
  integer               , intent(in)  :: nucleus_max_ang_mom(nucl_num)
  integer               , intent(in)  :: shell_ang_mom(shell_num)
  double precision      , intent(in)  :: ao_factor(ao_num)
  double precision      , intent(in)  :: shell_vgl(shell_num,5,point_num)
  double precision      , intent(out) :: ao_value(ao_num,point_num)

  double precision  :: e_coord(3), n_coord(3)
  integer*8         :: n_poly
  integer           :: l, il, k
  integer*8         :: ipoint, inucl, ishell
  integer*8         :: ishell_start, ishell_end
  integer           :: lstart(0:20)
  double precision  :: x, y, z, r2
  double precision  :: cutoff
  integer, external :: qmckl_ao_polynomial_vgl_doc_f

  double precision, allocatable  :: poly_vgl(:,:)
  integer         , allocatable  :: powers(:,:), ao_index(:)

  allocate(poly_vgl(5,ao_num), powers(3,ao_num), ao_index(ao_num))

  ! Pre-computed data
  do l=0,20
     lstart(l) = l*(l+1)*(l+2)/6 +1
  end do

  k=1
  do inucl=1,nucl_num
     ishell_start = nucleus_index(inucl) + 1
     ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
     do ishell = ishell_start, ishell_end
        l = shell_ang_mom(ishell)
        ao_index(ishell) = k
        k = k + lstart(l+1) - lstart(l)
     end do
  end do
  info = QMCKL_SUCCESS

  ! Don't compute polynomials when the radial part is zero.
  cutoff = 27.631021115928547  !-dlog(1.d-12)

  do ipoint = 1, point_num
     e_coord(1) = coord(ipoint,1)
     e_coord(2) = coord(ipoint,2)
     e_coord(3) = coord(ipoint,3)
     do inucl=1,nucl_num
        n_coord(1) = nucl_coord(inucl,1)
        n_coord(2) = nucl_coord(inucl,2)
        n_coord(3) = nucl_coord(inucl,3)

        ! Test if the point is in the range of the nucleus
        x = e_coord(1) - n_coord(1)
        y = e_coord(2) - n_coord(2)
        z = e_coord(3) - n_coord(3)

        r2 = x*x + y*y + z*z

        if (r2 > cutoff*nucleus_range(inucl)) then
           cycle
        end if

        ! Compute polynomials
        info = qmckl_ao_polynomial_vgl_doc_f(context, e_coord, n_coord, &
             nucleus_max_ang_mom(inucl), n_poly, powers, 3_8, &
             poly_vgl, 5_8)

        ! Loop over shells
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
        do ishell = ishell_start, ishell_end
           k = ao_index(ishell)
           l = shell_ang_mom(ishell)
           do il = lstart(l), lstart(l+1)-1
              ! Value
              ao_value(k,ipoint) = &
                   poly_vgl(1,il) * shell_vgl(ishell,1,ipoint) * ao_factor(k)
              k = k+1
           end do
        end do
     end do
  end do

  deallocate(poly_vgl, powers)
end function qmckl_compute_ao_value_doc_f



! #+CALL: generate_c_interface(table=qmckl_ao_value_args_doc,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_value_doc"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_ao_value_doc &
    (context, &
     ao_num, &
     shell_num, &
     point_num, &
     nucl_num, &
     coord, &
     nucl_coord, &
     nucleus_index, &
     nucleus_shell_num, &
     nucleus_range, &
     nucleus_max_ang_mom, &
     shell_ang_mom, &
     ao_factor, &
     shell_vgl, &
     ao_value) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int32_t) , intent(in)          :: nucleus_max_ang_mom(nucl_num)
  integer (c_int32_t) , intent(in)          :: shell_ang_mom(shell_num)
  real    (c_double ) , intent(in)          :: ao_factor(ao_num)
  real    (c_double ) , intent(in)          :: shell_vgl(shell_num,5,point_num)
  real    (c_double ) , intent(out)         :: ao_value(ao_num,point_num)

  integer(c_int32_t), external :: qmckl_compute_ao_value_doc_f
  info = qmckl_compute_ao_value_doc_f &
         (context, &
     ao_num, &
     shell_num, &
     point_num, &
     nucl_num, &
     coord, &
     nucl_coord, &
     nucleus_index, &
     nucleus_shell_num, &
     nucleus_range, &
     nucleus_max_ang_mom, &
     shell_ang_mom, &
     ao_factor, &
     shell_vgl, &
     ao_value)

end function qmckl_compute_ao_value_doc

! Unoptimized version
!      #+NAME: qmckl_ao_vgl_args_doc
!     | Variable              | Type                              | In/Out | Description                                  |
!     |-----------------------+-----------------------------------+--------+----------------------------------------------|
!     | ~context~             | ~qmckl_context~                   | in     | Global state                                 |
!     | ~ao_num~              | ~int64_t~                         | in     | Number of AOs                                |
!     | ~shell_num~           | ~int64_t~                         | in     | Number of shells                             |
!     | ~point_num~           | ~int64_t~                         | in     | Number of points                             |
!     | ~nucl_num~            | ~int64_t~                         | in     | Number of nuclei                             |
!     | ~coord~               | ~double[3][point_num]~            | in     | Coordinates                                  |
!     | ~nucl_coord~          | ~double[3][nucl_num]~             | in     | Nuclear  coordinates                         |
!     | ~nucleus_index~       | ~int64_t[nucl_num]~               | in     | Index of the 1st shell of each nucleus       |
!     | ~nucleus_shell_num~   | ~int64_t[nucl_num]~               | in     | Number of shells per nucleus                 |
!     | ~nucleus_range~       | ~double[nucl_num]~                | in     | Range beyond which all is zero               |
!     | ~nucleus_max_ang_mom~ | ~int32_t[nucl_num]~               | in     | Maximum angular momentum per nucleus         |
!     | ~shell_ang_mom~       | ~int32_t[shell_num]~              | in     | Angular momentum of each shell               |
!     | ~ao_factor~           | ~double[ao_num]~                  | in     | Normalization factor of the AOs              |
!     | ~shell_vgl~           | ~double[point_num][5][shell_num]~ | in     | Value, gradients and Laplacian of the shells |
!     | ~ao_vgl~              | ~double[point_num][5][ao_num]~    | out    | Value, gradients and Laplacian of the AOs    |


integer function qmckl_compute_ao_vgl_doc_f(context, &
     ao_num, shell_num, point_num, nucl_num, &
     coord, nucl_coord, nucleus_index, nucleus_shell_num, &
     nucleus_range, nucleus_max_ang_mom, shell_ang_mom, &
     ao_factor, shell_vgl, ao_vgl) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: ao_num
  integer*8             , intent(in)  :: shell_num
  integer*8             , intent(in)  :: point_num
  integer*8             , intent(in)  :: nucl_num
  double precision      , intent(in)  :: coord(point_num,3)
  double precision      , intent(in)  :: nucl_coord(nucl_num,3)
  integer*8             , intent(in)  :: nucleus_index(nucl_num)
  integer*8             , intent(in)  :: nucleus_shell_num(nucl_num)
  double precision      , intent(in)  :: nucleus_range(nucl_num)
  integer               , intent(in)  :: nucleus_max_ang_mom(nucl_num)
  integer               , intent(in)  :: shell_ang_mom(shell_num)
  double precision      , intent(in)  :: ao_factor(ao_num)
  double precision      , intent(in)  :: shell_vgl(shell_num,5,point_num)
  double precision      , intent(out) :: ao_vgl(ao_num,5,point_num)

  double precision  :: e_coord(3), n_coord(3)
  integer*8         :: n_poly
  integer           :: l, il, k
  integer*8         :: ipoint, inucl, ishell
  integer*8         :: ishell_start, ishell_end
  integer           :: lstart(0:20)
  double precision  :: x, y, z, r2
  double precision  :: cutoff
  integer, external :: qmckl_ao_polynomial_vgl_doc_f

  double precision, allocatable  :: poly_vgl(:,:)
  integer         , allocatable  :: powers(:,:), ao_index(:)

  allocate(poly_vgl(5,ao_num), powers(3,ao_num), ao_index(ao_num))

  ! Pre-computed data
  do l=0,20
     lstart(l) = l*(l+1)*(l+2)/6 +1
  end do

  k=1
  do inucl=1,nucl_num
     ishell_start = nucleus_index(inucl) + 1
     ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
     do ishell = ishell_start, ishell_end
        l = shell_ang_mom(ishell)
        ao_index(ishell) = k
        k = k + lstart(l+1) - lstart(l)
     end do
  end do
  info = QMCKL_SUCCESS

  ! Don't compute polynomials when the radial part is zero.
  cutoff = 27.631021115928547  ! -dlog(1.d-12)

  do ipoint = 1, point_num
     e_coord(1) = coord(ipoint,1)
     e_coord(2) = coord(ipoint,2)
     e_coord(3) = coord(ipoint,3)
     do inucl=1,nucl_num
        n_coord(1) = nucl_coord(inucl,1)
        n_coord(2) = nucl_coord(inucl,2)
        n_coord(3) = nucl_coord(inucl,3)

        ! Test if the point is in the range of the nucleus
        x = e_coord(1) - n_coord(1)
        y = e_coord(2) - n_coord(2)
        z = e_coord(3) - n_coord(3)

        r2 = x*x + y*y + z*z

        if (r2 > cutoff*nucleus_range(inucl)) then
           cycle
        end if

        ! Compute polynomials
        info = qmckl_ao_polynomial_vgl_doc_f(context, e_coord, n_coord, &
             nucleus_max_ang_mom(inucl), n_poly, powers, 3_8, &
             poly_vgl, 5_8)

        ! Loop over shells
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
        do ishell = ishell_start, ishell_end
           k = ao_index(ishell)
           l = shell_ang_mom(ishell)
           do il = lstart(l), lstart(l+1)-1
              ! Value
              ao_vgl(k,1,ipoint) = &
                   poly_vgl(1,il) * shell_vgl(ishell,1,ipoint) * ao_factor(k)

              ! Grad_x
              ao_vgl(k,2,ipoint) = ( &
                   poly_vgl(2,il) * shell_vgl(ishell,1,ipoint) + &
                   poly_vgl(1,il) * shell_vgl(ishell,2,ipoint) &
                   ) * ao_factor(k)

              ! Grad_y
              ao_vgl(k,3,ipoint) = ( &
                   poly_vgl(3,il) * shell_vgl(ishell,1,ipoint) + &
                   poly_vgl(1,il) * shell_vgl(ishell,3,ipoint) &
                   ) * ao_factor(k)

              ! Grad_z
              ao_vgl(k,4,ipoint) = ( &
                   poly_vgl(4,il) * shell_vgl(ishell,1,ipoint) + &
                   poly_vgl(1,il) * shell_vgl(ishell,4,ipoint) &
                   ) * ao_factor(k)

              ! Lapl_z
              ao_vgl(k,5,ipoint) = ( &
                   poly_vgl(5,il) * shell_vgl(ishell,1,ipoint) + &
                   poly_vgl(1,il) * shell_vgl(ishell,5,ipoint) + &
                   2.d0 * ( &
                   poly_vgl(2,il) * shell_vgl(ishell,2,ipoint) + &
                   poly_vgl(3,il) * shell_vgl(ishell,3,ipoint) + &
                   poly_vgl(4,il) * shell_vgl(ishell,4,ipoint) ) &
                   ) * ao_factor(k)

              k = k+1
           end do
        end do
     end do
  end do

  deallocate(poly_vgl, powers)
end function qmckl_compute_ao_vgl_doc_f



! #+CALL: generate_c_interface(table=qmckl_ao_vgl_args_doc,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_vgl_doc"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_ao_vgl_doc &
    (context, &
     ao_num, &
     shell_num, &
     point_num, &
     nucl_num, &
     coord, &
     nucl_coord, &
     nucleus_index, &
     nucleus_shell_num, &
     nucleus_range, &
     nucleus_max_ang_mom, &
     shell_ang_mom, &
     ao_factor, &
     shell_vgl, &
     ao_vgl) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int32_t) , intent(in)          :: nucleus_max_ang_mom(nucl_num)
  integer (c_int32_t) , intent(in)          :: shell_ang_mom(shell_num)
  real    (c_double ) , intent(in)          :: ao_factor(ao_num)
  real    (c_double ) , intent(in)          :: shell_vgl(shell_num,5,point_num)
  real    (c_double ) , intent(out)         :: ao_vgl(ao_num,5,point_num)

  integer(c_int32_t), external :: qmckl_compute_ao_vgl_doc_f
  info = qmckl_compute_ao_vgl_doc_f &
         (context, &
     ao_num, &
     shell_num, &
     point_num, &
     nucl_num, &
     coord, &
     nucl_coord, &
     nucleus_index, &
     nucleus_shell_num, &
     nucleus_range, &
     nucleus_max_ang_mom, &
     shell_ang_mom, &
     ao_factor, &
     shell_vgl, &
     ao_vgl)

end function qmckl_compute_ao_vgl_doc
