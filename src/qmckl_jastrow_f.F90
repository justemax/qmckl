! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_asymp_jasb
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_asymp_jasb_args
!     | Variable                  | Type                 | In/Out | Description             |
!     |---------------------------+----------------------+--------+-------------------------|
!     | ~context~                 | ~qmckl_context~      | in     | Global state            |
!     | ~bord_num~                | ~int64_t~            | in     | Order of the polynomial |
!     | ~bord_vector~             | ~double[bord_num+1]~ | in     | Values of b             |
!     | ~rescale_factor_kappa_ee~ | ~double~             | in     | Electron coordinates    |
!     | ~asymp_jasb~              | ~double[2]~          | out    | Asymptotic value        |


integer function qmckl_compute_asymp_jasb_f(context, bord_num, bord_vector, rescale_factor_kappa_ee, asymp_jasb) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: bord_num
  double precision      , intent(in)  :: bord_vector(bord_num + 1)
  double precision      , intent(in)  :: rescale_factor_kappa_ee
  double precision      , intent(out) :: asymp_jasb(2)

  integer*8 :: i, p
  double precision   :: kappa_inv, x, asym_one
  kappa_inv = 1.0d0 / rescale_factor_kappa_ee

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (bord_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  asym_one = bord_vector(1) * kappa_inv / (1.0d0 + bord_vector(2) * kappa_inv)
  asymp_jasb(:) = (/asym_one, 0.5d0 * asym_one/)

  do i = 1, 2
     x = kappa_inv
     do p = 2, bord_num
        x = x * kappa_inv
        asymp_jasb(i) = asymp_jasb(i) + bord_vector(p + 1) * x
     end do
  end do

end function qmckl_compute_asymp_jasb_f

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_ee
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_ee_args
!     | Variable               | Type                                   | In/Out | Description                 |
!     |------------------------+----------------------------------------+--------+-----------------------------|
!     | ~context~              | ~qmckl_context~                        | in     | Global state                |
!     | ~walk_num~             | ~int64_t~                              | in     | Number of walkers           |
!     | ~elec_num~             | ~int64_t~                              | in     | Number of electrons         |
!     | ~up_num~               | ~int64_t~                              | in     | Number of alpha electrons   |
!     | ~bord_num~             | ~int64_t~                              | in     | Number of coefficients      |
!     | ~bord_vector~          | ~double[bord_num+1]~                   | in     | List of coefficients        |
!     | ~ee_distance_rescaled~ | ~double[walk_num][elec_num][elec_num]~ | in     | Electron-electron distances |
!     | ~asymp_jasb~           | ~double[2]~                            | in     | Electron-electron distances |
!     | ~factor_ee~            | ~double[walk_num]~                     | out    | Electron-electron distances |


integer function qmckl_compute_factor_ee_f(context, walk_num, elec_num, up_num, bord_num,            &
                                           bord_vector, ee_distance_rescaled, asymp_jasb, factor_ee) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num, elec_num, bord_num, up_num
  double precision      , intent(in)  :: bord_vector(bord_num + 1)
  double precision      , intent(in)  :: ee_distance_rescaled(elec_num, elec_num, walk_num)
  double precision      , intent(in)  :: asymp_jasb(2)
  double precision      , intent(out) :: factor_ee(walk_num)

  integer*8 :: i, j, p, ipar, nw
  double precision   :: x, power_ser, spin_fact

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (bord_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  factor_ee = 0.0d0

  do nw =1, walk_num
  do j = 1, elec_num
     do i = 1, j - 1
        x = ee_distance_rescaled(i,j,nw)
        power_ser = 0.0d0
        spin_fact = 1.0d0
        ipar = 1

        do p = 2, bord_num
          x = x * ee_distance_rescaled(i,j,nw)
          power_ser = power_ser + bord_vector(p + 1) * x
        end do

        if(j <= up_num .OR. i > up_num) then
          spin_fact = 0.5d0
          ipar = 2
        endif

        factor_ee(nw) = factor_ee(nw) + spin_fact * bord_vector(1)  * &
                                ee_distance_rescaled(i,j,nw) / &
                                (1.0d0 + bord_vector(2) *   &
                                ee_distance_rescaled(i,j,nw))  &
                               -asymp_jasb(ipar) + power_ser

     end do
  end do
  end do

end function qmckl_compute_factor_ee_f

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_ee_deriv_e
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_ee_deriv_e_args
!     | Variable                       | Type                                      | In/Out | Description                 |
!     |--------------------------------+-------------------------------------------+--------+-----------------------------|
!     | ~context~                      | ~qmckl_context~                           | in     | Global state                |
!     | ~walk_num~                     | ~int64_t~                                 | in     | Number of walkers           |
!     | ~elec_num~                     | ~int64_t~                                 | in     | Number of electrons         |
!     | ~up_num~                       | ~int64_t~                                 | in     | Number of alpha electrons   |
!     | ~bord_num~                     | ~int64_t~                                 | in     | Number of coefficients      |
!     | ~bord_vector~                  | ~double[bord_num+1]~                      | in     | List of coefficients        |
!     | ~ee_distance_rescaled~         | ~double[walk_num][elec_num][elec_num]~    | in     | Electron-electron distances |
!     | ~ee_distance_rescaled_deriv_e~ | ~double[walk_num][4][elec_num][elec_num]~ | in     | Electron-electron distances |
!     | ~factor_ee_deriv_e~            | ~double[walk_num][4][elec_num]~           | out    | Electron-electron distances |


integer function qmckl_compute_factor_ee_deriv_e_doc_f( &
     context, walk_num, elec_num, up_num, bord_num, &
     bord_vector, ee_distance_rescaled, ee_distance_rescaled_deriv_e,  &
     factor_ee_deriv_e) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num, elec_num, bord_num, up_num
  double precision      , intent(in)  :: bord_vector(bord_num + 1)
  double precision      , intent(in)  :: ee_distance_rescaled(elec_num, elec_num,walk_num)
  double precision      , intent(in)  :: ee_distance_rescaled_deriv_e(4,elec_num, elec_num,walk_num)   !TODO
  double precision      , intent(out) :: factor_ee_deriv_e(elec_num,4,walk_num)

  integer*8 :: i, j, p, nw, ii
  double precision   :: x, spin_fact, y
  double precision   :: den, invden, invden2, invden3, xinv
  double precision   :: lap1, lap2, lap3, third
  double precision, dimension(3) :: pow_ser_g
  double precision, dimension(4) :: dx

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (bord_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  factor_ee_deriv_e = 0.0d0
  third = 1.0d0 / 3.0d0

  do nw =1, walk_num
  do j = 1, elec_num
     do i = 1, elec_num
        x = ee_distance_rescaled(i,j,nw)
        if(abs(x) < 1.0d-18) cycle
        pow_ser_g   = 0.0d0
        spin_fact   = 1.0d0
        den         = 1.0d0 + bord_vector(2) * x
        invden      = 1.0d0 / den
        invden2     = invden * invden
        invden3     = invden2 * invden
        xinv        = 1.0d0 / (x + 1.0d-18)

        dx(1) = ee_distance_rescaled_deriv_e(1, i, j, nw)
        dx(2) = ee_distance_rescaled_deriv_e(2, i, j, nw)
        dx(3) = ee_distance_rescaled_deriv_e(3, i, j, nw)
        dx(4) = ee_distance_rescaled_deriv_e(4, i, j, nw)

        if((i .LE. up_num .AND. j .LE. up_num ) .OR.  &
           (i .GT. up_num .AND. j .GT. up_num)) then
          spin_fact = 0.5d0
        endif

        lap1 = 0.0d0
        lap2 = 0.0d0
        lap3 = 0.0d0
        do ii = 1, 3
          x = ee_distance_rescaled(i, j, nw)
          if(abs(x) < 1.0d-18) cycle
          do p = 2, bord_num
            y = p * bord_vector(p + 1) * x
            pow_ser_g(ii) = pow_ser_g(ii) + y * dx(ii)
            lap1 = lap1 + (p - 1) * y * xinv * dx(ii) * dx(ii)
            lap2 = lap2 + y
            x = x * ee_distance_rescaled(i, j, nw)
          end do

          lap3 = lap3 - 2.0d0 * bord_vector(2) * dx(ii) * dx(ii)

          factor_ee_deriv_e( j, ii, nw) = factor_ee_deriv_e( j, ii, nw) + spin_fact * bord_vector(1)  * &
                                dx(ii) * invden2 + pow_ser_g(ii)
        end do

        ii = 4
        lap2 = lap2 * dx(ii) * third
        lap3 = lap3 + den * dx(ii)
        lap3 = lap3 * (spin_fact * bord_vector(1) * invden3)
        factor_ee_deriv_e( j, ii, nw) = factor_ee_deriv_e( j, ii, nw) + lap1 + lap2 + lap3

     end do
  end do
  end do

end function qmckl_compute_factor_ee_deriv_e_doc_f




! #+CALL: generate_c_interface(table=qmckl_factor_ee_deriv_e_args,rettyp=get_value("CRetType"),fname="qmckl_compute_factor_ee_deriv_e_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_factor_ee_deriv_e_doc &
    (context, &
         walk_num, &
         elec_num, &
         up_num, &
         bord_num, &
         bord_vector, &
         ee_distance_rescaled, &
         ee_distance_rescaled_deriv_e, &
         factor_ee_deriv_e) &
        bind(C) result(info)

      use, intrinsic :: iso_c_binding
      implicit none

      integer (c_int64_t) , intent(in)  , value :: context
      integer (c_int64_t) , intent(in)  , value :: walk_num
      integer (c_int64_t) , intent(in)  , value :: elec_num
      integer (c_int64_t) , intent(in)  , value :: up_num
      integer (c_int64_t) , intent(in)  , value :: bord_num
      real    (c_double ) , intent(in)          :: bord_vector(bord_num+1)
      real    (c_double ) , intent(in)          :: ee_distance_rescaled(elec_num,elec_num,walk_num)
      real    (c_double ) , intent(in)          :: ee_distance_rescaled_deriv_e(elec_num,elec_num,4,walk_num)
      real    (c_double ) , intent(out)         :: factor_ee_deriv_e(elec_num,4,walk_num)

      integer(c_int32_t), external :: qmckl_compute_factor_ee_deriv_e_doc_f
      info = qmckl_compute_factor_ee_deriv_e_doc_f &
         (context, &
         walk_num, &
         elec_num, &
         up_num, &
         bord_num, &
         bord_vector, &
         ee_distance_rescaled, &
         ee_distance_rescaled_deriv_e, &
         factor_ee_deriv_e)

    end function qmckl_compute_factor_ee_deriv_e_doc

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_en
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_en_args
!     | Variable               | Type                                   | In/Out | Description                |
!     |------------------------+----------------------------------------+--------+----------------------------|
!     | ~context~              | ~qmckl_context~                        | in     | Global state               |
!     | ~walk_num~             | ~int64_t~                              | in     | Number of walkers          |
!     | ~elec_num~             | ~int64_t~                              | in     | Number of electrons        |
!     | ~nucl_num~             | ~int64_t~                              | in     | Number of nucleii          |
!     | ~type_nucl_num~        | ~int64_t~                              | in     | Number of unique nuclei    |
!     | ~type_nucl_vector~     | ~int64_t[nucl_num]~                    | in     | IDs of unique nucleii      |
!     | ~aord_num~             | ~int64_t~                              | in     | Number of coefficients     |
!     | ~aord_vector~          | ~double[aord_num+1][type_nucl_num]~    | in     | List of coefficients       |
!     | ~en_distance_rescaled~ | ~double[walk_num][nucl_num][elec_num]~ | in     | Electron-nucleus distances |
!     | ~factor_en~            | ~double[walk_num]~                     | out    | Electron-nucleus jastrow   |


integer function qmckl_compute_factor_en_f( &
     context, walk_num, elec_num, nucl_num, type_nucl_num, &
     type_nucl_vector, aord_num, aord_vector, &
     en_distance_rescaled, factor_en) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num, elec_num, aord_num, nucl_num, type_nucl_num
  integer*8             , intent(in)  :: type_nucl_vector(nucl_num)
  double precision      , intent(in)  :: aord_vector(aord_num + 1, type_nucl_num)
  double precision      , intent(in)  :: en_distance_rescaled(elec_num, nucl_num, walk_num)
  double precision      , intent(out) :: factor_en(walk_num)

  integer*8 :: i, a, p, nw
  double precision   :: x, power_ser

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (aord_num <= 0) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  factor_en = 0.0d0

  do nw =1, walk_num
  do a = 1, nucl_num
     do i = 1, elec_num
        x = en_distance_rescaled(i, a, nw)
        power_ser = 0.0d0

        do p = 2, aord_num
          x = x * en_distance_rescaled(i, a, nw)
          power_ser = power_ser + aord_vector(p + 1, type_nucl_vector(a)) * x
        end do

        factor_en(nw) = factor_en(nw) + aord_vector(1, type_nucl_vector(a)) *  &
                                en_distance_rescaled(i, a, nw) /               &
                                (1.0d0 + aord_vector(2, type_nucl_vector(a)) * &
                                en_distance_rescaled(i, a, nw))                &
                                + power_ser

     end do
  end do
  end do

end function qmckl_compute_factor_en_f

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_en_deriv_e
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_en_deriv_e_args
!     | Variable                       | Type                                      | In/Out | Description                           |
!     |--------------------------------+-------------------------------------------+--------+---------------------------------------|
!     | ~context~                      | ~qmckl_context~                           | in     | Global state                          |
!     | ~walk_num~                     | ~int64_t~                                 | in     | Number of walkers                     |
!     | ~elec_num~                     | ~int64_t~                                 | in     | Number of electrons                   |
!     | ~nucl_num~                     | ~int64_t~                                 | in     | Number of nucleii                     |
!     | ~type_nucl_num~                | ~int64_t~                                 | in     | Number of unique nuclei               |
!     | ~type_nucl_vector~             | ~int64_t[nucl_num]~                       | in     | IDs of unique nucleii                 |
!     | ~aord_num~                     | ~int64_t~                                 | in     | Number of coefficients                |
!     | ~aord_vector~                  | ~double[aord_num+1][type_nucl_num]~       | in     | List of coefficients                  |
!     | ~en_distance_rescaled~         | ~double[walk_num][nucl_num][elec_num]~    | in     | Electron-nucleus distances            |
!     | ~en_distance_rescaled_deriv_e~ | ~double[walk_num][4][nucl_num][elec_num]~ | in     | Electron-nucleus distance derivatives |
!     | ~factor_en_deriv_e~            | ~double[walk_num][4][elec_num]~           | out    | Electron-nucleus jastrow              |


integer function qmckl_compute_factor_en_deriv_e_f( &
     context, walk_num, elec_num, nucl_num, type_nucl_num, &
     type_nucl_vector, aord_num, aord_vector, &
     en_distance_rescaled, en_distance_rescaled_deriv_e, factor_en_deriv_e) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num, elec_num, aord_num, nucl_num, type_nucl_num
  integer*8             , intent(in)  :: type_nucl_vector(nucl_num)
  double precision      , intent(in)  :: aord_vector(aord_num + 1, type_nucl_num)
  double precision      , intent(in)  :: en_distance_rescaled(elec_num, nucl_num, walk_num)
  double precision      , intent(in)  :: en_distance_rescaled_deriv_e(4, elec_num, nucl_num, walk_num)
  double precision      , intent(out) :: factor_en_deriv_e(elec_num,4,walk_num)

  integer*8 :: i, a, p, ipar, nw, ii
  double precision   :: x, den, invden, invden2, invden3, xinv
  double precision   :: y, lap1, lap2, lap3, third
  double precision, dimension(3) :: power_ser_g
  double precision, dimension(4) :: dx

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (aord_num <= 0) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  factor_en_deriv_e = 0.0d0
  third = 1.0d0 / 3.0d0

  do nw =1, walk_num
  do a = 1, nucl_num
     do i = 1, elec_num
        x = en_distance_rescaled(i,a,nw)
        if(abs(x) < 1.0d-18) continue
        power_ser_g = 0.0d0
        den = 1.0d0 + aord_vector(2, type_nucl_vector(a)) * x
        invden = 1.0d0 / den
        invden2 = invden * invden
        invden3 = invden2 * invden
        xinv = 1.0d0 / x

        do ii = 1, 4
          dx(ii) = en_distance_rescaled_deriv_e(ii,i,a,nw)
        end do

        lap1 = 0.0d0
        lap2 = 0.0d0
        lap3 = 0.0d0
        do ii = 1, 3
          x = en_distance_rescaled(i, a, nw)
          do p = 2, aord_num
            y = p * aord_vector(p + 1, type_nucl_vector(a)) * x
            power_ser_g(ii) = power_ser_g(ii) + y * dx(ii)
            lap1 = lap1 + (p - 1) * y * xinv * dx(ii) * dx(ii)
            lap2 = lap2 + y
            x = x * en_distance_rescaled(i, a, nw)
          end do

          lap3 = lap3 - 2.0d0 * aord_vector(2, type_nucl_vector(a)) * dx(ii) * dx(ii)

          factor_en_deriv_e(i, ii, nw) = factor_en_deriv_e(i, ii, nw) + aord_vector(1, type_nucl_vector(a)) &
                                  * dx(ii) * invden2                                                        &
                                  + power_ser_g(ii)

        end do

        ii = 4
        lap2 = lap2 * dx(ii) * third
        lap3 = lap3 + den * dx(ii)
        lap3 = lap3 * aord_vector(1, type_nucl_vector(a)) * invden3
        factor_en_deriv_e(i, ii, nw) = factor_en_deriv_e(i, ii, nw) + lap1 + lap2 + lap3

     end do
  end do
  end do

end function qmckl_compute_factor_en_deriv_e_f




! #+CALL: generate_c_interface(table=qmckl_factor_en_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_factor_en_deriv_e &
    (context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     type_nucl_num, &
     type_nucl_vector, &
     aord_num, &
     aord_vector, &
     en_distance_rescaled, &
     en_distance_rescaled_deriv_e, &
     factor_en_deriv_e) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  integer (c_int64_t) , intent(in)  , value :: aord_num
  real    (c_double ) , intent(in)          :: aord_vector(type_nucl_num,aord_num+1)
  real    (c_double ) , intent(in)          :: en_distance_rescaled(elec_num,nucl_num,walk_num)
  real    (c_double ) , intent(in)          :: en_distance_rescaled_deriv_e(elec_num,nucl_num,4,walk_num)
  real    (c_double ) , intent(out)         :: factor_en_deriv_e(elec_num,4,walk_num)

  integer(c_int32_t), external :: qmckl_compute_factor_en_deriv_e_f
  info = qmckl_compute_factor_en_deriv_e_f &
         (context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     type_nucl_num, &
     type_nucl_vector, &
     aord_num, &
     aord_vector, &
     en_distance_rescaled, &
     en_distance_rescaled_deriv_e, &
     factor_en_deriv_e)

end function qmckl_compute_factor_en_deriv_e

! CPU

!    :PROPERTIES:
!    :Name:     qmckl_compute_een_rescaled_e
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_een_rescaled_e_args
!     | Variable                  | Type                                               | In/Out | Description                          |
!     |---------------------------+----------------------------------------------------+--------+--------------------------------------|
!     | ~context~                 | ~qmckl_context~                                    | in     | Global state                         |
!     | ~walk_num~                | ~int64_t~                                          | in     | Number of walkers                    |
!     | ~elec_num~                | ~int64_t~                                          | in     | Number of electrons                  |
!     | ~cord_num~                | ~int64_t~                                          | in     | Order of polynomials                 |
!     | ~rescale_factor_kappa_ee~ | ~double~                                           | in     | Factor to rescale ee distances       |
!     | ~ee_distance~             | ~double[walk_num][elec_num][elec_num]~             | in     | Electron-electron distances          |
!     | ~een_rescaled_e~          | ~double[walk_num][0:cord_num][elec_num][elec_num]~ | out    | Electron-electron rescaled distances |


integer function qmckl_compute_een_rescaled_e_doc_f( &
     context, walk_num, elec_num, cord_num, rescale_factor_kappa_ee,  &
     ee_distance, een_rescaled_e) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: cord_num
  double precision      , intent(in)  :: rescale_factor_kappa_ee
  double precision      , intent(in)  :: ee_distance(elec_num,elec_num,walk_num)
  double precision      , intent(out) :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  double precision,dimension(:,:),allocatable :: een_rescaled_e_ij
  double precision                    :: x
  integer*8                           :: i, j, k, l, nw

  allocate(een_rescaled_e_ij(elec_num * (elec_num - 1) / 2, cord_num + 1))

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  ! Prepare table of exponentiated distances raised to appropriate power
  een_rescaled_e             = 0.0d0
  do nw = 1, walk_num
  een_rescaled_e_ij       = 0.0d0
  een_rescaled_e_ij(:, 1) = 1.0d0


  k = 0
  do j = 1, elec_num
    do i = 1, j - 1
      k = k + 1
      een_rescaled_e_ij(k, 2) = dexp(-rescale_factor_kappa_ee * ee_distance(i, j, nw))
    end do
  end do


  do l = 2, cord_num
    do k = 1, elec_num * (elec_num - 1)/2
      een_rescaled_e_ij(k, l + 1) = een_rescaled_e_ij(k, l + 1 - 1) * een_rescaled_e_ij(k, 2)
    end do
  end do

  ! prepare the actual een table
  een_rescaled_e(:, :, 0, nw) = 1.0d0

  do l = 1, cord_num
    k = 0
    do j = 1, elec_num
      do i = 1, j - 1
        k = k + 1
        x = een_rescaled_e_ij(k, l + 1)
        een_rescaled_e(i, j, l, nw) = x
        een_rescaled_e(j, i, l, nw) = x
      end do
    end do
  end do

  do l = 0, cord_num
    do j = 1, elec_num
      een_rescaled_e(j, j, l, nw) = 0.0d0
    end do
  end do

  end do

end function qmckl_compute_een_rescaled_e_doc_f



! #+CALL: generate_c_interface(table=qmckl_factor_een_rescaled_e_args,rettyp=get_value("CRetType"),fname="qmckl_compute_een_rescaled_e_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_een_rescaled_e_doc &
    (context, walk_num, elec_num, cord_num, rescale_factor_kappa_ee, &
    ee_distance, een_rescaled_e) &
bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_kappa_ee
  real    (c_double ) , intent(in)          :: ee_distance(elec_num,elec_num,walk_num)
  real    (c_double ) , intent(out)         :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_een_rescaled_e_doc_f
  info = qmckl_compute_een_rescaled_e_doc_f &
	 (context, walk_num, elec_num, cord_num, rescale_factor_kappa_ee, ee_distance, een_rescaled_e)

end function qmckl_compute_een_rescaled_e_doc

! CPU

!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_een_rescaled_e_deriv_e
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_een_rescaled_e_deriv_e_args
!     | Variable                  | Type                                                  | In/Out | Description                          |
!     |---------------------------+-------------------------------------------------------+--------+--------------------------------------|
!     | ~context~                 | ~qmckl_context~                                       | in     | Global state                         |
!     | ~walk_num~                | ~int64_t~                                             | in     | Number of walkers                    |
!     | ~elec_num~                | ~int64_t~                                             | in     | Number of electrons                  |
!     | ~cord_num~                | ~int64_t~                                             | in     | Order of polynomials                 |
!     | ~rescale_factor_kappa_ee~ | ~double~                                              | in     | Factor to rescale ee distances       |
!     | ~coord_ee~                | ~double[walk_num][3][elec_num]~                       | in     | Electron coordinates                 |
!     | ~ee_distance~             | ~double[walk_num][elec_num][elec_num]~                | in     | Electron-electron distances          |
!     | ~een_rescaled_e~          | ~double[walk_num][0:cord_num][elec_num][elec_num]~    | in     | Electron-electron distances          |
!     | ~een_rescaled_e_deriv_e~  | ~double[walk_num][0:cord_num][elec_num][4][elec_num]~ | out    | Electron-electron rescaled distances |


integer function qmckl_compute_factor_een_rescaled_e_deriv_e_f( &
     context, walk_num, elec_num, cord_num, rescale_factor_kappa_ee,  &
     coord_ee, ee_distance, een_rescaled_e, een_rescaled_e_deriv_e) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: cord_num
  double precision      , intent(in)  :: rescale_factor_kappa_ee
  double precision      , intent(in)  :: coord_ee(elec_num,3,walk_num)
  double precision      , intent(in)  :: ee_distance(elec_num,elec_num,walk_num)
  double precision      , intent(in)  :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  double precision      , intent(out) :: een_rescaled_e_deriv_e(elec_num,4,elec_num,0:cord_num,walk_num)
  double precision,dimension(:,:,:),allocatable  :: elec_dist_deriv_e
  double precision                    :: x, rij_inv, kappa_l
  integer*8                           :: i, j, k, l, nw, ii

  allocate(elec_dist_deriv_e(4,elec_num,elec_num))

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  ! Prepare table of exponentiated distances raised to appropriate power
  een_rescaled_e_deriv_e     = 0.0d0
  do nw = 1, walk_num
    do j = 1, elec_num
      do i = 1, elec_num
        rij_inv = 1.0d0 / ee_distance(i, j, nw)
        do ii = 1, 3
          elec_dist_deriv_e(ii, i, j) = (coord_ee(i, ii, nw) - coord_ee(j, ii, nw)) * rij_inv
        end do
        elec_dist_deriv_e(4, i, j) = 2.0d0 * rij_inv
      end do
      elec_dist_deriv_e(:, j, j) = 0.0d0
    end do

    ! prepare the actual een table
    do l = 1, cord_num
      kappa_l = - dble(l) * rescale_factor_kappa_ee
      do j = 1, elec_num
        do i = 1, elec_num
          een_rescaled_e_deriv_e(i, 1, j, l, nw) = kappa_l * elec_dist_deriv_e(1, i, j)
          een_rescaled_e_deriv_e(i, 2, j, l, nw) = kappa_l * elec_dist_deriv_e(2, i, j)
          een_rescaled_e_deriv_e(i, 3, j, l, nw) = kappa_l * elec_dist_deriv_e(3, i, j)
          een_rescaled_e_deriv_e(i, 4, j, l, nw) = kappa_l * elec_dist_deriv_e(4, i, j)

          een_rescaled_e_deriv_e(i, 4, j, l, nw) = een_rescaled_e_deriv_e(i, 4, j, l, nw)              &
                    + een_rescaled_e_deriv_e(i, 1, j, l, nw) * een_rescaled_e_deriv_e(i, 1, j, l, nw)  &
                    + een_rescaled_e_deriv_e(i, 2, j, l, nw) * een_rescaled_e_deriv_e(i, 2, j, l, nw)  &
                    + een_rescaled_e_deriv_e(i, 3, j, l, nw) * een_rescaled_e_deriv_e(i, 3, j, l, nw)

          een_rescaled_e_deriv_e(i, 1, j, l, nw) = een_rescaled_e_deriv_e(i, 1, j, l, nw) *   &
                                                    een_rescaled_e(i, j, l, nw)
          een_rescaled_e_deriv_e(i, 3, j, l, nw) = een_rescaled_e_deriv_e(i, 2, j, l, nw) *   &
                                                    een_rescaled_e(i, j, l, nw)
          een_rescaled_e_deriv_e(i, 3, j, l, nw) = een_rescaled_e_deriv_e(i, 3, j, l, nw) *   &
                                                    een_rescaled_e(i, j, l, nw)
          een_rescaled_e_deriv_e(i, 4, j, l, nw) = een_rescaled_e_deriv_e(i, 4, j, l, nw) *   &
                                                    een_rescaled_e(i, j, l, nw)
        end do
      end do
    end do
  end do

end function qmckl_compute_factor_een_rescaled_e_deriv_e_f




! #+CALL: generate_c_interface(table=qmckl_factor_een_rescaled_e_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_factor_een_rescaled_e_deriv_e &
    (context, &
     walk_num, &
     elec_num, &
     cord_num, &
     rescale_factor_kappa_ee, &
     coord_ee, &
     ee_distance, &
     een_rescaled_e, &
     een_rescaled_e_deriv_e) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_kappa_ee
  real    (c_double ) , intent(in)          :: coord_ee(elec_num,3,walk_num)
  real    (c_double ) , intent(in)          :: ee_distance(elec_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: een_rescaled_e_deriv_e(elec_num,4,elec_num,0:cord_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_factor_een_rescaled_e_deriv_e_f
  info = qmckl_compute_factor_een_rescaled_e_deriv_e_f &
         (context, &
     walk_num, &
     elec_num, &
     cord_num, &
     rescale_factor_kappa_ee, &
     coord_ee, &
     ee_distance, &
     een_rescaled_e, &
     een_rescaled_e_deriv_e)

end function qmckl_compute_factor_een_rescaled_e_deriv_e

! CPU

!    :PROPERTIES:
!    :Name:     qmckl_compute_een_rescaled_n
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_een_rescaled_n_args
!     | Variable                  | Type                                               | In/Out | Description                         |
!     |---------------------------+----------------------------------------------------+--------+-------------------------------------|
!     | ~context~                 | ~qmckl_context~                                    | in     | Global state                        |
!     | ~walk_num~                | ~int64_t~                                          | in     | Number of walkers                   |
!     | ~elec_num~                | ~int64_t~                                          | in     | Number of electrons                 |
!     | ~nucl_num~                | ~int64_t~                                          | in     | Number of atoms                     |
!     | ~cord_num~                | ~int64_t~                                          | in     | Order of polynomials                |
!     | ~rescale_factor_kappa_en~ | ~double~                                           | in     | Factor to rescale ee distances      |
!     | ~en_distance~             | ~double[walk_num][elec_num][nucl_num]~             | in     | Electron-nucleus distances          |
!     | ~een_rescaled_n~          | ~double[walk_num][0:cord_num][nucl_num][elec_num]~ | out    | Electron-nucleus rescaled distances |


integer function qmckl_compute_een_rescaled_n_f( &
     context, walk_num, elec_num, nucl_num, cord_num, rescale_factor_kappa_en,  &
     en_distance, een_rescaled_n) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: nucl_num
  integer*8             , intent(in)  :: cord_num
  double precision      , intent(in)  :: rescale_factor_kappa_en
  double precision      , intent(in)  :: en_distance(elec_num,nucl_num,walk_num)
  double precision      , intent(out) :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  double precision                    :: x
  integer*8                           :: i, a, k, l, nw

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  ! Prepare table of exponentiated distances raised to appropriate power
  een_rescaled_n             = 0.0d0
  do nw = 1, walk_num

  ! prepare the actual een table
  een_rescaled_n(:, :, 0, nw) = 1.0d0

  do a = 1, nucl_num
    do i = 1, elec_num
      een_rescaled_n(i, a, 1, nw) = dexp(-rescale_factor_kappa_en * en_distance(i, a, nw))
    end do
  end do

  do l = 2, cord_num
    do a = 1, nucl_num
      do i = 1, elec_num
        een_rescaled_n(i, a, l, nw) = een_rescaled_n(i, a, l - 1, nw) * een_rescaled_n(i, a, 1, nw)
      end do
    end do
  end do
  end do

end function qmckl_compute_een_rescaled_n_f

! CPU

!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_een_rescaled_n_deriv_e
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_compute_factor_een_rescaled_n_deriv_e_args
!     | Variable                  | Type                                                  | In/Out | Description                         |
!     |---------------------------+-------------------------------------------------------+--------+-------------------------------------|
!     | ~context~                 | ~qmckl_context~                                       | in     | Global state                        |
!     | ~walk_num~                | ~int64_t~                                             | in     | Number of walkers                   |
!     | ~elec_num~                | ~int64_t~                                             | in     | Number of electrons                 |
!     | ~nucl_num~                | ~int64_t~                                             | in     | Number of atoms                     |
!     | ~cord_num~                | ~int64_t~                                             | in     | Order of polynomials                |
!     | ~rescale_factor_kappa_en~ | ~double~                                              | in     | Factor to rescale ee distances      |
!     | ~coord_ee~                | ~double[walk_num][3][elec_num]~                       | in     | Electron coordinates                |
!     | ~coord_en~                | ~double[3][nucl_num]~                                 | in     | Nuclear coordinates                 |
!     | ~en_distance~             | ~double[walk_num][elec_num][nucl_num]~                | in     | Electron-nucleus distances          |
!     | ~een_rescaled_n~          | ~double[walk_num][0:cord_num][nucl_num][elec_num]~    | in     | Electron-nucleus distances          |
!     | ~een_rescaled_n_deriv_e~  | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~ | out    | Electron-nucleus rescaled distances |


integer function qmckl_compute_factor_een_rescaled_n_deriv_e_f( &
     context, walk_num, elec_num, nucl_num, &
     cord_num, rescale_factor_kappa_en, &
     coord_ee, coord_en, en_distance, een_rescaled_n, een_rescaled_n_deriv_e) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: nucl_num
  integer*8             , intent(in)  :: cord_num
  double precision      , intent(in)  :: rescale_factor_kappa_en
  double precision      , intent(in)  :: coord_ee(elec_num,3,walk_num)
  double precision      , intent(in)  :: coord_en(nucl_num,3)
  double precision      , intent(in)  :: en_distance(elec_num,nucl_num,walk_num)
  double precision      , intent(in)  :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  double precision      , intent(out) :: een_rescaled_n_deriv_e(elec_num,4,nucl_num,0:cord_num,walk_num)
  double precision,dimension(:,:,:),allocatable :: elnuc_dist_deriv_e
  double precision                    :: x, ria_inv, kappa_l
  integer*8                           :: i, a, k, l, nw, ii

  allocate(elnuc_dist_deriv_e(4, elec_num, nucl_num))

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  ! Prepare table of exponentiated distances raised to appropriate power
  een_rescaled_n_deriv_e             = 0.0d0
  do nw = 1, walk_num

  ! prepare the actual een table
  do a = 1, nucl_num
    do i = 1, elec_num
      ria_inv = 1.0d0 / en_distance(i, a, nw)
      do ii = 1, 3
        elnuc_dist_deriv_e(ii, i, a) = (coord_ee(i, ii, nw) - coord_en(a, ii)) * ria_inv
      end do
      elnuc_dist_deriv_e(4, i, a) = 2.0d0 * ria_inv
    end do
  end do

  do l = 0, cord_num
    kappa_l = - dble(l) * rescale_factor_kappa_en
    do a = 1, nucl_num
      do i = 1, elec_num
        een_rescaled_n_deriv_e(i, 1, a, l, nw) = kappa_l * elnuc_dist_deriv_e(1, i, a)
        een_rescaled_n_deriv_e(i, 2, a, l, nw) = kappa_l * elnuc_dist_deriv_e(2, i, a)
        een_rescaled_n_deriv_e(i, 3, a, l, nw) = kappa_l * elnuc_dist_deriv_e(3, i, a)
        een_rescaled_n_deriv_e(i, 4, a, l, nw) = kappa_l * elnuc_dist_deriv_e(4, i, a)

        een_rescaled_n_deriv_e(i, 4, a, l, nw) = een_rescaled_n_deriv_e(i, 4, a, l, nw)           &
                + een_rescaled_n_deriv_e(i, 1, a, l, nw) * een_rescaled_n_deriv_e(i, 1, a, l, nw) &
                + een_rescaled_n_deriv_e(i, 2, a, l, nw) * een_rescaled_n_deriv_e(i, 2, a, l, nw) &
                + een_rescaled_n_deriv_e(i, 3, a, l, nw) * een_rescaled_n_deriv_e(i, 3, a, l, nw)

        een_rescaled_n_deriv_e(i, 1, a, l, nw) = een_rescaled_n_deriv_e(i, 1, a, l, nw) * &
                                                  een_rescaled_n(i, a, l, nw)
        een_rescaled_n_deriv_e(i, 2, a, l, nw) = een_rescaled_n_deriv_e(i, 2, a, l, nw) * &
                                                  een_rescaled_n(i, a, l, nw)
        een_rescaled_n_deriv_e(i, 3, a, l, nw) = een_rescaled_n_deriv_e(i, 3, a, l, nw) * &
                                                  een_rescaled_n(i, a, l, nw)
        een_rescaled_n_deriv_e(i, 4, a, l, nw) = een_rescaled_n_deriv_e(i, 4, a, l, nw) * &
                                                  een_rescaled_n(i, a, l, nw)
      end do
    end do
  end do
  end do

end function qmckl_compute_factor_een_rescaled_n_deriv_e_f



! #+CALL: generate_c_interface(table=qmckl_compute_factor_een_rescaled_n_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_factor_een_rescaled_n_deriv_e &
    (context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     rescale_factor_kappa_en, &
     coord_ee, &
     coord_en, &
     en_distance, &
     een_rescaled_n, &
     een_rescaled_n_deriv_e) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_kappa_en
  real    (c_double ) , intent(in)          :: coord_ee(elec_num,3,walk_num)
  real    (c_double ) , intent(in)          :: coord_en(nucl_num,3)
  real    (c_double ) , intent(in)          :: en_distance(nucl_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: een_rescaled_n_deriv_e(elec_num,4,nucl_num,0:cord_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_factor_een_rescaled_n_deriv_e_f
  info = qmckl_compute_factor_een_rescaled_n_deriv_e_f &
         (context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     rescale_factor_kappa_en, &
     coord_ee, &
     coord_en, &
     en_distance, &
     een_rescaled_n, &
     een_rescaled_n_deriv_e)

end function qmckl_compute_factor_een_rescaled_n_deriv_e

! Compute dim_cord_vect

!    :PROPERTIES:
!    :Name:     qmckl_compute_dim_cord_vect
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_dim_cord_vect_args
!     | Variable        | Type            | In/Out | Description                       |
!     |-----------------+-----------------+--------+-----------------------------------|
!     | ~context~       | ~qmckl_context~ | in     | Global state                      |
!     | ~cord_num~      | ~int64_t~       | in     | Order of polynomials              |
!     | ~dim_cord_vect~ | ~int64_t~       | out    | dimension of cord_vect_full table |


integer function qmckl_compute_dim_cord_vect_f( &
     context, cord_num, dim_cord_vect) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: cord_num
  integer*8             , intent(out) :: dim_cord_vect
  double precision                    :: x
  integer*8                           :: i, a, k, l, p, lmax

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  dim_cord_vect = 0

  do p = 2, cord_num
    do k = p - 1, 0, -1
      if (k .ne. 0) then
        lmax = p - k
      else
        lmax = p - k - 2
      endif
      do l = lmax, 0, -1
        if (iand(p - k - l, 1_8) == 1) cycle
        dim_cord_vect = dim_cord_vect + 1
      end do
    end do
  end do

end function qmckl_compute_dim_cord_vect_f

! CPU & offload version
!    :PROPERTIES:
!    :Name:     qmckl_compute_cord_vect_full
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_cord_vect_full_args
!     | Variable           | Type                                   | In/Out | Description                  |
!     |--------------------+----------------------------------------+--------+------------------------------|
!     | ~context~          | ~qmckl_context~                        | in     | Global state                 |
!     | ~nucl_num~         | ~int64_t~                              | in     | Number of atoms              |
!     | ~dim_cord_vect~    | ~int64_t~                              | in     | dimension of cord full table |
!     | ~type_nucl_num~    | ~int64_t~                              | in     | dimension of cord full table |
!     | ~type_nucl_vector~ | ~int64_t[nucl_num]~                    | in     | dimension of cord full table |
!     | ~cord_vector~      | ~double[dim_cord_vect][type_nucl_num]~ | in     | dimension of cord full table |
!     | ~cord_vect_full~   | ~double[dim_cord_vect][nucl_num]~      | out    | Full list of coefficients    |


integer function qmckl_compute_cord_vect_full_doc_f( &
     context, nucl_num, dim_cord_vect, type_nucl_num,  &
     type_nucl_vector, cord_vector, cord_vect_full) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: nucl_num
  integer*8             , intent(in)  :: dim_cord_vect
  integer*8             , intent(in)  :: type_nucl_num
  integer*8             , intent(in)  :: type_nucl_vector(nucl_num)
  double precision      , intent(in)  :: cord_vector(type_nucl_num, dim_cord_vect)
  double precision      , intent(out) :: cord_vect_full(nucl_num,dim_cord_vect)
  double precision                    :: x
  integer*8                           :: i, a, k, l, nw

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (type_nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (dim_cord_vect <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif


  do a = 1, nucl_num
    cord_vect_full(a,1:dim_cord_vect) = cord_vector(type_nucl_vector(a),1:dim_cord_vect)
  end do

end function qmckl_compute_cord_vect_full_doc_f



! #+CALL: generate_c_interface(table=qmckl_factor_cord_vect_full_args,rettyp=get_value("CRetType"),fname="qmckl_compute_cord_vect_full_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_cord_vect_full_doc &
    (context, nucl_num, dim_cord_vect, type_nucl_num, type_nucl_vector, cord_vector, cord_vect_full) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: dim_cord_vect
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  real    (c_double ) , intent(in)          :: cord_vector(type_nucl_num,dim_cord_vect)
  real    (c_double ) , intent(out)         :: cord_vect_full(nucl_num,dim_cord_vect)

  integer(c_int32_t), external :: qmckl_compute_cord_vect_full_doc_f
  info = qmckl_compute_cord_vect_full_doc_f &
	 (context, nucl_num, dim_cord_vect, type_nucl_num, type_nucl_vector, cord_vector, cord_vect_full)

end function qmckl_compute_cord_vect_full_doc

! CPU & offload version
!    :PROPERTIES:
!    :Name:     qmckl_compute_lkpm_combined_index
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_lkpm_combined_index_args
!     | Variable              | Type                        | In/Out | Description                   |
!     |-----------------------+-----------------------------+--------+-------------------------------|
!     | ~context~             | ~qmckl_context~             | in     | Global state                  |
!     | ~cord_num~            | ~int64_t~                   | in     | Order of polynomials          |
!     | ~dim_cord_vect~       | ~int64_t~                   | in     | dimension of cord full table  |
!     | ~lkpm_combined_index~ | ~int64_t[4][dim_cord_vect]~ | out    | Full list of combined indices |


integer function qmckl_compute_lkpm_combined_index_f( &
     context, cord_num, dim_cord_vect,  lkpm_combined_index) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: cord_num
  integer*8             , intent(in)  :: dim_cord_vect
  integer*8             , intent(out) :: lkpm_combined_index(dim_cord_vect, 4)
  double precision                    :: x
  integer*8                           :: i, a, k, l, kk, p, lmax, m

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (dim_cord_vect <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif


  kk = 0
  do p = 2, cord_num
    do k = p - 1, 0, -1
      if (k .ne. 0) then
        lmax = p - k
      else
        lmax = p - k - 2
      end if
      do l = lmax, 0, -1
        if (iand(p - k - l, 1_8) .eq. 1) cycle
        m = (p - k - l)/2
        kk = kk + 1
        lkpm_combined_index(kk, 1) = l
        lkpm_combined_index(kk, 2) = k
        lkpm_combined_index(kk, 3) = p
        lkpm_combined_index(kk, 4) = m
      end do
    end do
  end do

end function qmckl_compute_lkpm_combined_index_f

integer function qmckl_compute_tmp_c_doc_f( &
     context, cord_num, elec_num, nucl_num, &
     walk_num, een_rescaled_e, een_rescaled_n, tmp_c) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: cord_num
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: nucl_num
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: een_rescaled_e(elec_num, elec_num, 0:cord_num, walk_num)
  double precision      , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  double precision      , intent(out) :: tmp_c(elec_num, nucl_num,0:cord_num, 0:cord_num-1, walk_num)
  double precision                    :: x
  integer*8                           :: i, j, a, l, kk, p, lmax, nw
  character                           :: TransA, TransB
  double precision                    :: alpha, beta
  integer*8                           :: M, N, K, LDA, LDB, LDC

  TransA = 'N'
  TransB = 'N'
  alpha = 1.0d0
  beta  = 0.0d0

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  M = elec_num
  N = nucl_num*(cord_num + 1)
  K = elec_num
  LDA = size(een_rescaled_e,1)
  LDB = size(een_rescaled_n,1)
  LDC = size(tmp_c,1)

  do nw=1, walk_num
  do i=0, cord_num-1
  info = qmckl_dgemm(context, TransA, TransB, M, N, K, alpha,     &
                     een_rescaled_e(1,1,i,nw),LDA*1_8,                     &
                     een_rescaled_n(1,1,0,nw),LDB*1_8,                     &
                     beta,                                       &
                     tmp_c(1,1,0,i,nw),LDC)
  end do
  end do

end function qmckl_compute_tmp_c_doc_f



!     #+CALL: generate_c_interface(table=qmckl_factor_tmp_c_args,rettyp=get_value("FRetType"),fname="qmckl_compute_tmp_c_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_tmp_c_doc &
    (context, cord_num, elec_num, nucl_num, walk_num, een_rescaled_e, een_rescaled_n, tmp_c) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: tmp_c(elec_num,nucl_num,0:cord_num,0:cord_num-1,walk_num)

  integer(c_int32_t), external :: qmckl_compute_tmp_c_doc_f
  info = qmckl_compute_tmp_c_doc_f &
         (context, cord_num, elec_num, nucl_num, walk_num, een_rescaled_e, een_rescaled_n, tmp_c)

end function qmckl_compute_tmp_c_doc

integer function qmckl_compute_dtmp_c_doc_f( &
     context, cord_num, elec_num, nucl_num, &
     walk_num, een_rescaled_e_deriv_e, een_rescaled_n, dtmp_c) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: cord_num
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: nucl_num
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: een_rescaled_e_deriv_e(elec_num, 4, elec_num, 0:cord_num, walk_num)
  double precision      , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  double precision      , intent(out) :: dtmp_c(elec_num, 4, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  double precision                    :: x
  integer*8                           :: i, j, a, l, kk, p, lmax, nw, ii
  character                           :: TransA, TransB
  double precision                    :: alpha, beta
  integer*8                           :: M, N, K, LDA, LDB, LDC

  TransA = 'N'
  TransB = 'N'
  alpha = 1.0d0
  beta  = 0.0d0

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  M = 4*elec_num
  N = nucl_num*(cord_num + 1)
  K = elec_num
  LDA = 4*size(een_rescaled_e_deriv_e,1)
  LDB = size(een_rescaled_n,1)
  LDC = 4*size(dtmp_c,1)

  do nw=1, walk_num
     do i=0, cord_num-1
        info = qmckl_dgemm(context,TransA, TransB, M, N, K, alpha,  &
             een_rescaled_e_deriv_e(1,1,1,i,nw),LDA*1_8,            &
             een_rescaled_n(1,1,0,nw),LDB*1_8,                      &
             beta,                                                  &
             dtmp_c(1,1,1,0,i,nw),LDC)
     end do
  end do

end function qmckl_compute_dtmp_c_doc_f



!     #+CALL: generate_c_interface(table=qmckl_factor_dtmp_c_args,rettyp=get_value("FRetType"),fname="qmckl_compute_dtmp_c_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_dtmp_c_doc &
    (context, cord_num, elec_num, nucl_num, walk_num, een_rescaled_e_deriv_e, een_rescaled_n, dtmp_c) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: een_rescaled_e_deriv_e(elec_num,4,elec_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: dtmp_c(elec_num,nucl_num,0:cord_num,0:cord_num-1,walk_num)

  integer(c_int32_t), external :: qmckl_compute_dtmp_c_doc_f
  info = qmckl_compute_dtmp_c_doc_f &
         (context, cord_num, elec_num, nucl_num, walk_num, een_rescaled_e_deriv_e, een_rescaled_n, dtmp_c)

end function qmckl_compute_dtmp_c_doc

! Compute naive
!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_een_naive
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_een_naive_args
!     | Variable              | Type                                               | In/Out | Description                          |
!     |-----------------------+----------------------------------------------------+--------+--------------------------------------|
!     | ~context~             | ~qmckl_context~                                    | in     | Global state                         |
!     | ~walk_num~            | ~int64_t~                                          | in     | Number of walkers                    |
!     | ~elec_num~            | ~int64_t~                                          | in     | Number of electrons                  |
!     | ~nucl_num~            | ~int64_t~                                          | in     | Number of nucleii                    |
!     | ~cord_num~            | ~int64_t~                                          | in     | order of polynomials                 |
!     | ~dim_cord_vect~       | ~int64_t~                                          | in     | dimension of full coefficient vector |
!     | ~cord_vect_full~      | ~double[dim_cord_vect][nucl_num]~                  | in     | full coefficient vector              |
!     | ~lkpm_combined_index~ | ~int64_t[4][dim_cord_vect]~                        | in     | combined indices                     |
!     | ~een_rescaled_e~      | ~double[walk_num][elec_num][elec_num][0:cord_num]~ | in     | Electron-nucleus rescaled            |
!     | ~een_rescaled_n~      | ~double[walk_num][elec_num][nucl_num][0:cord_num]~ | in     | Electron-nucleus rescaled factor     |
!     | ~factor_een~          | ~double[walk_num]~                                 | out    | Electron-nucleus jastrow             |


integer function qmckl_compute_factor_een_naive_f( &
     context, walk_num, elec_num, nucl_num, cord_num,&
     dim_cord_vect, cord_vect_full, lkpm_combined_index, &
     een_rescaled_e, een_rescaled_n, factor_een) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num, elec_num, cord_num, nucl_num, dim_cord_vect
  integer*8             , intent(in)  :: lkpm_combined_index(dim_cord_vect,4)
  double precision      , intent(in)  :: cord_vect_full(nucl_num, dim_cord_vect)
  double precision      , intent(in)  :: een_rescaled_e(0:cord_num, elec_num, elec_num, walk_num)
  double precision      , intent(in)  :: een_rescaled_n(0:cord_num, nucl_num, elec_num, walk_num)
  double precision      , intent(out) :: factor_een(walk_num)

  integer*8 :: i, a, j, l, k, p, m, n, nw
  double precision :: accu, accu2, cn

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  factor_een = 0.0d0

  do nw =1, walk_num
  do n = 1, dim_cord_vect
    l = lkpm_combined_index(n, 1)
    k = lkpm_combined_index(n, 2)
    p = lkpm_combined_index(n, 3)
    m = lkpm_combined_index(n, 4)

    do a = 1, nucl_num
      accu2 = 0.0d0
      cn = cord_vect_full(a, n)
      do j = 1, elec_num
        accu = 0.0d0
        do i = 1, elec_num
          accu = accu + een_rescaled_e(k,i,j,nw) *       &
                        een_rescaled_n(m,a,i,nw)
          !if(nw .eq. 1) then
          !  print *,l,k,p,m,j,i,een_rescaled_e(k,i,j,nw), een_rescaled_n(m,a,i,nw), accu
          !endif
        end do
        accu2 = accu2 + accu * een_rescaled_n(m + l,a,j,nw)
        !print *, l,m,nw,accu, accu2, een_rescaled_n(m + l, a, j, nw), cn, factor_een(nw)
      end do
      factor_een(nw) = factor_een(nw) + accu2 * cn
    end do
  end do
  end do

end function qmckl_compute_factor_een_naive_f



! #+CALL: generate_c_interface(table=qmckl_factor_een_naive_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_factor_een_naive &
(context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     dim_cord_vect, &
     cord_vect_full, &
     lkpm_combined_index, &
     een_rescaled_e, &
     een_rescaled_n, &
     factor_een) &
bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_cord_vect
  real    (c_double ) , intent(in)          :: cord_vect_full(nucl_num,dim_cord_vect)
  integer (c_int64_t) , intent(in)          :: lkpm_combined_index(dim_cord_vect,4)
  real    (c_double ) , intent(in)          :: een_rescaled_e(0:cord_num,elec_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(0:cord_num,nucl_num,elec_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_een(walk_num)

  integer(c_int32_t), external :: qmckl_compute_factor_een_naive_f
  info = qmckl_compute_factor_een_naive_f &
     (context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     dim_cord_vect, &
     cord_vect_full, &
     lkpm_combined_index, &
     een_rescaled_e, &
     een_rescaled_n, &
     factor_een)

end function qmckl_compute_factor_een_naive

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_een
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_een_args
!     | Variable              | Type                                                             | In/Out                          | Description                          |
!     |-----------------------+------------------------------------------------------------------+---------------------------------+--------------------------------------|
!     | ~context~             | ~qmckl_context~                                                  | in                              | Global state                         |
!     | ~walk_num~            | ~int64_t~                                                        | in                              | Number of walkers                    |
!     | ~elec_num~            | ~int64_t~                                                        | in                              | Number of electrons                  |
!     | ~nucl_num~            | ~int64_t~                                                        | in                              | Number of nucleii                    |
!     | ~cord_num~            | ~int64_t~                                                        | in                              | order of polynomials                 |
!     | ~dim_cord_vect~       | ~int64_t~                                                        | in                              | dimension of full coefficient vector |
!     | ~cord_vect_full~      | ~double[dim_cord_vect][nucl_num]~                                | in                              | full coefficient vector              |
!     | ~lkpm_combined_index~ | ~int64_t[4][dim_cord_vect]~                                      | in                              | combined indices                     |
!     | ~tmp_c~               | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~ | vector of non-zero coefficients |                                      |
!     | ~een_rescaled_n~      | ~double[walk_num][0:cord_num][nucl_num][elec_num]~               | in                              | Electron-nucleus rescaled factor     |
!     | ~factor_een~          | ~double[walk_num]~                                               | out                             | Electron-nucleus jastrow             |


integer function qmckl_compute_factor_een_f( &
     context, walk_num, elec_num, nucl_num, cord_num,   &
     dim_cord_vect, cord_vect_full, lkpm_combined_index, &
     tmp_c, een_rescaled_n, factor_een) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num, elec_num, cord_num, nucl_num, dim_cord_vect
  integer*8             , intent(in)  :: lkpm_combined_index(dim_cord_vect,4)
  double precision      , intent(in)  :: cord_vect_full(nucl_num, dim_cord_vect)
  double precision      , intent(in)  :: tmp_c(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  double precision      , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  double precision      , intent(out) :: factor_een(walk_num)

  integer*8 :: i, a, j, l, k, p, m, n, nw
  double precision :: accu, accu2, cn

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  factor_een = 0.0d0

  do nw =1, walk_num
  do n = 1, dim_cord_vect
    l = lkpm_combined_index(n, 1)
    k = lkpm_combined_index(n, 2)
    p = lkpm_combined_index(n, 3)
    m = lkpm_combined_index(n, 4)

    do a = 1, nucl_num
      cn = cord_vect_full(a, n)
      if(cn == 0.d0) cycle

      accu = 0.0d0
      do j = 1, elec_num
        accu = accu + een_rescaled_n(j,a,m,nw) * tmp_c(j,a,m+l,k,nw)
      end do
      factor_een(nw) = factor_een(nw) + accu * cn
    end do
  end do
  end do

end function qmckl_compute_factor_een_f



! #+CALL: generate_c_interface(table=qmckl_factor_een_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_factor_een &
(context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     dim_cord_vect, &
     cord_vect_full, &
     lkpm_combined_index, &
     een_rescaled_e, &
     een_rescaled_n, &
     factor_een) &
bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_cord_vect
  real    (c_double ) , intent(in)          :: cord_vect_full(nucl_num,dim_cord_vect)
  integer (c_int64_t) , intent(in)          :: lkpm_combined_index(dim_cord_vect,4)
  real    (c_double ) , intent(in)          :: een_rescaled_e(0:cord_num,elec_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_een(walk_num)

  integer(c_int32_t), external :: qmckl_compute_factor_een_f
  info = qmckl_compute_factor_een_f &
     (context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     dim_cord_vect, &
     cord_vect_full, &
     lkpm_combined_index, &
     een_rescaled_e, &
     een_rescaled_n, &
     factor_een)

end function qmckl_compute_factor_een

! Compute Naive
!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_een_deriv_e_naive
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_een_deriv_e_naive_args
!     | Variable                 | Type                                                  | In/Out | Description                          |
!     |--------------------------+-------------------------------------------------------+--------+--------------------------------------|
!     | ~context~                | ~qmckl_context~                                       | in     | Global state                         |
!     | ~walk_num~               | ~int64_t~                                             | in     | Number of walkers                    |
!     | ~elec_num~               | ~int64_t~                                             | in     | Number of electrons                  |
!     | ~nucl_num~               | ~int64_t~                                             | in     | Number of nucleii                    |
!     | ~cord_num~               | ~int64_t~                                             | in     | order of polynomials                 |
!     | ~dim_cord_vect~          | ~int64_t~                                             | in     | dimension of full coefficient vector |
!     | ~cord_vect_full~         | ~double[dim_cord_vect][nucl_num]~                     | in     | full coefficient vector              |
!     | ~lkpm_combined_index~    | ~int64_t[4][dim_cord_vect]~                           | in     | combined indices                     |
!     | ~een_rescaled_e~         | ~double[walk_num][elec_num][elec_num][0:cord_num]~    | in     | Electron-nucleus rescaled            |
!     | ~een_rescaled_n~         | ~double[walk_num][elec_num][nucl_num][0:cord_num]~    | in     | Electron-nucleus rescaled factor     |
!     | ~een_rescaled_e_deriv_e~ | ~double[walk_num][elec_num][4][elec_num][0:cord_num]~ | in     | Electron-nucleus rescaled            |
!     | ~een_rescaled_n_deriv_e~ | ~double[walk_num][elec_num][4][nucl_num][0:cord_num]~ | in     | Electron-nucleus rescaled factor     |
!     | ~factor_een_deriv_e~     | ~double[walk_num][4][elec_num]~                       | out    | Electron-nucleus jastrow             |


integer function qmckl_compute_factor_een_deriv_e_naive_f( &
     context, walk_num, elec_num, nucl_num, cord_num, dim_cord_vect, &
     cord_vect_full, lkpm_combined_index, een_rescaled_e, een_rescaled_n, &
     een_rescaled_e_deriv_e, een_rescaled_n_deriv_e, factor_een_deriv_e)&
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num, elec_num, cord_num, nucl_num, dim_cord_vect
  integer*8             , intent(in)  :: lkpm_combined_index(dim_cord_vect, 4)
  double precision      , intent(in)  :: cord_vect_full(nucl_num, dim_cord_vect)
  double precision      , intent(in)  :: een_rescaled_e(0:cord_num, elec_num, elec_num, walk_num)
  double precision      , intent(in)  :: een_rescaled_n(0:cord_num, nucl_num, elec_num, walk_num)
  double precision      , intent(in)  :: een_rescaled_e_deriv_e(0:cord_num, elec_num, 4, elec_num, walk_num)
  double precision      , intent(in)  :: een_rescaled_n_deriv_e(0:cord_num, nucl_num, 4, elec_num, walk_num)
  double precision      , intent(out) :: factor_een_deriv_e(elec_num, 4, walk_num)

  integer*8 :: i, a, j, l, k, p, m, n, nw
  double precision :: accu, accu2, cn
  double precision :: daccu(1:4), daccu2(1:4)

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  factor_een_deriv_e = 0.0d0

  do nw =1, walk_num
  do n = 1, dim_cord_vect
    l = lkpm_combined_index(n, 1)
    k = lkpm_combined_index(n, 2)
    p = lkpm_combined_index(n, 3)
    m = lkpm_combined_index(n, 4)

    do a = 1, nucl_num
      cn = cord_vect_full(a, n)
      do j = 1, elec_num
        accu = 0.0d0
        accu2 = 0.0d0
        daccu = 0.0d0
        daccu2 = 0.0d0
        do i = 1, elec_num
          accu = accu + een_rescaled_e(k, i, j, nw) *         &
                        een_rescaled_n(m, a, i, nw)
          accu2 = accu2 + een_rescaled_e(k, i, j, nw) *       &
                          een_rescaled_n(m + l, a, i, nw)
          daccu(1:4) = daccu(1:4) + een_rescaled_e_deriv_e(k, j, 1:4, i, nw) *         &
                                    een_rescaled_n(m, a, i, nw)
          daccu2(1:4) = daccu2(1:4) + een_rescaled_e_deriv_e(k, j, 1:4, i, nw) *       &
                                      een_rescaled_n(m + l, a, i, nw)
        end do
        factor_een_deriv_e(j, 1:4, nw) = factor_een_deriv_e(j, 1:4, nw) +              &
               (accu * een_rescaled_n_deriv_e(m + l, a, 1:4, j, nw)                    &
                + daccu(1:4) * een_rescaled_n(m + l, a, j, nw)                         &
                + daccu2(1:4) * een_rescaled_n(m, a, j, nw)                            &
                + accu2 * een_rescaled_n_deriv_e(m, a, 1:4, j, nw)) * cn

        factor_een_deriv_e(j, 4, nw) = factor_een_deriv_e(j, 4, nw) + 2.0d0 * (        &
            daccu (1) * een_rescaled_n_deriv_e(m + l, a, 1, j, nw) +                    &
            daccu (2) * een_rescaled_n_deriv_e(m + l, a, 2, j, nw) +                    &
            daccu (3) * een_rescaled_n_deriv_e(m + l, a, 3, j, nw) +                    &
            daccu2(1) * een_rescaled_n_deriv_e(m, a, 1, j, nw    ) +                    &
            daccu2(2) * een_rescaled_n_deriv_e(m, a, 2, j, nw    ) +                    &
            daccu2(3) * een_rescaled_n_deriv_e(m, a, 3, j, nw    ) ) * cn

      end do
    end do
  end do
  end do

end function qmckl_compute_factor_een_deriv_e_naive_f




! #+CALL: generate_c_interface(table=qmckl_factor_een_deriv_e_naive_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_factor_een_deriv_e_naive &
(context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     dim_cord_vect, &
     cord_vect_full, &
     lkpm_combined_index, &
     een_rescaled_e, &
     een_rescaled_n, &
     een_rescaled_e_deriv_e, &
     een_rescaled_n_deriv_e, &
     factor_een_deriv_e) &
bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_cord_vect
  real    (c_double ) , intent(in)          :: cord_vect_full(nucl_num,dim_cord_vect)
  integer (c_int64_t) , intent(in)          :: lkpm_combined_index(dim_cord_vect,4)
  real    (c_double ) , intent(in)          :: een_rescaled_e(0:cord_num,elec_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(0:cord_num,nucl_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_e_deriv_e(0:cord_num,elec_num,4,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n_deriv_e(0:cord_num,nucl_num,4,elec_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_een_deriv_e(elec_num,4,walk_num)

  integer(c_int32_t), external :: qmckl_compute_factor_een_deriv_e_naive_f
  info = qmckl_compute_factor_een_deriv_e_naive_f &
     (context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     dim_cord_vect, &
     cord_vect_full, &
     lkpm_combined_index, &
     een_rescaled_e, &
     een_rescaled_n, &
     een_rescaled_e_deriv_e, &
     een_rescaled_n_deriv_e, &
     factor_een_deriv_e)

end function qmckl_compute_factor_een_deriv_e_naive

! CPU

!    :PROPERTIES:
!    :Name:     qmckl_compute_factor_een_deriv_e
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_factor_een_deriv_e_args
!     | Variable                 | Type                                                                | In/Out | Description                                    |
!     |--------------------------+---------------------------------------------------------------------+--------+------------------------------------------------|
!     | ~context~                | ~qmckl_context~                                                     | in     | Global state                                   |
!     | ~walk_num~               | ~int64_t~                                                           | in     | Number of walkers                              |
!     | ~elec_num~               | ~int64_t~                                                           | in     | Number of electrons                            |
!     | ~nucl_num~               | ~int64_t~                                                           | in     | Number of nucleii                              |
!     | ~cord_num~               | ~int64_t~                                                           | in     | order of polynomials                           |
!     | ~dim_cord_vect~          | ~int64_t~                                                           | in     | dimension of full coefficient vector           |
!     | ~cord_vect_full~         | ~double[dim_cord_vect][nucl_num]~                                   | in     | full coefficient vector                        |
!     | ~lkpm_combined_index~    | ~int64_t[4][dim_cord_vect]~                                         | in     | combined indices                               |
!     | ~tmp_c~                  | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~    | in     | Temporary intermediate tensor                  |
!     | ~dtmp_c~                 | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][4][elec_num]~ | in     | vector of non-zero coefficients                |
!     | ~een_rescaled_n~         | ~double[walk_num][0:cord_num][nucl_num][elec_num]~                  | in     | Electron-nucleus rescaled factor               |
!     | ~een_rescaled_n_deriv_e~ | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~               | in     | Derivative of Electron-nucleus rescaled factor |
!     | ~factor_een_deriv_e~     | ~double[walk_num][4][elec_num]~                                     | out    | Derivative of Electron-nucleus jastrow         |



integer function qmckl_compute_factor_een_deriv_e_f( &
     context, walk_num, elec_num, nucl_num, &
     cord_num, dim_cord_vect, cord_vect_full, lkpm_combined_index, &
     tmp_c, dtmp_c, een_rescaled_n, een_rescaled_n_deriv_e, factor_een_deriv_e)&
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: walk_num, elec_num, cord_num, nucl_num, dim_cord_vect
  integer*8             , intent(in)  :: lkpm_combined_index(dim_cord_vect,4)
  double precision      , intent(in)  :: cord_vect_full(nucl_num, dim_cord_vect)
  double precision      , intent(in)  :: tmp_c(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  double precision      , intent(in)  :: dtmp_c(elec_num, 4, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  double precision      , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  double precision      , intent(in)  :: een_rescaled_n_deriv_e(elec_num, 4, nucl_num, 0:cord_num, walk_num)
  double precision      , intent(out) :: factor_een_deriv_e(elec_num,4,walk_num)

  integer*8 :: i, a, j, l, k, p, m, n, nw, ii
  double precision :: accu, accu2, cn

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (cord_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  factor_een_deriv_e = 0.0d0

  do nw =1, walk_num
  do n = 1, dim_cord_vect
    l = lkpm_combined_index(n, 1)
    k = lkpm_combined_index(n, 2)
    p = lkpm_combined_index(n, 3)
    m = lkpm_combined_index(n, 4)

    do a = 1, nucl_num
      cn = cord_vect_full(a, n)
      if(cn == 0.d0) cycle

      do ii = 1, 4
        do j = 1, elec_num
          factor_een_deriv_e(j,ii,nw) = factor_een_deriv_e(j,ii,nw) +                           (&
                                        tmp_c(j,a,m,k,nw)       * een_rescaled_n_deriv_e(j,ii,a,m+l,nw) + &
                                        (dtmp_c(j,ii,a,m,k,nw))   * een_rescaled_n(j,a,m+l,nw)            + &
                                        (dtmp_c(j,ii,a,m+l,k,nw)) * een_rescaled_n(j,a,m  ,nw)              + &
                                        tmp_c(j,a,m+l,k,nw)     * een_rescaled_n_deriv_e(j,ii,a,m,nw)     &
                                        ) * cn
        end do
      end do

      cn = cn + cn
      do j = 1, elec_num
        factor_een_deriv_e(j,4,nw) = factor_een_deriv_e(j,4,nw) +                              (&
                                      (dtmp_c(j,1,a,m  ,k,nw)) * een_rescaled_n_deriv_e(j,1,a,m+l,nw)  + &
                                      (dtmp_c(j,2,a,m  ,k,nw)) * een_rescaled_n_deriv_e(j,2,a,m+l,nw)  + &
                                      (dtmp_c(j,3,a,m  ,k,nw)) * een_rescaled_n_deriv_e(j,3,a,m+l,nw)  + &
                                      (dtmp_c(j,1,a,m+l,k,nw)) * een_rescaled_n_deriv_e(j,1,a,m  ,nw)  + &
                                      (dtmp_c(j,2,a,m+l,k,nw)) * een_rescaled_n_deriv_e(j,2,a,m  ,nw)  + &
                                      (dtmp_c(j,3,a,m+l,k,nw)) * een_rescaled_n_deriv_e(j,3,a,m  ,nw)    &
                                      ) * cn
      end do
    end do
  end do
  end do

end function qmckl_compute_factor_een_deriv_e_f




! #+CALL: generate_c_interface(table=qmckl_factor_een_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_factor_een_deriv_e &
(context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     dim_cord_vect, &
     cord_vect_full, &
     lkpm_combined_index, &
     tmp_c, &
     dtmp_c, &
     een_rescaled_n, &
     een_rescaled_n_deriv_e, &
     factor_een_deriv_e) &
bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_cord_vect
  real    (c_double ) , intent(in)          :: cord_vect_full(nucl_num,dim_cord_vect)
  integer (c_int64_t) , intent(in)          :: lkpm_combined_index(dim_cord_vect,4)
  real    (c_double ) , intent(in)          :: tmp_c(elec_num,nucl_num,0:cord_num,0:cord_num-1,walk_num)
  real    (c_double ) , intent(in)          :: dtmp_c(elec_num,4,nucl_num,0:cord_num,0:cord_num-1,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n_deriv_e(elec_num,4,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_een_deriv_e(elec_num,4,walk_num)

  integer(c_int32_t), external :: qmckl_compute_factor_een_deriv_e_f
  info = qmckl_compute_factor_een_deriv_e_f &
     (context, &
     walk_num, &
     elec_num, &
     nucl_num, &
     cord_num, &
     dim_cord_vect, &
     cord_vect_full, &
     lkpm_combined_index, &
     tmp_c, &
     dtmp_c, &
     een_rescaled_n, &
     een_rescaled_n_deriv_e, &
     factor_een_deriv_e)

end function qmckl_compute_factor_een_deriv_e
