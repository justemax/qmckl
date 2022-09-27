! CPU

!    :PROPERTIES:
!    :Name:     qmckl_compute_ee_distance
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_ee_distance_args
!     | Variable      | Type                                   | In/Out | Description                 |
!     |---------------+----------------------------------------+--------+-----------------------------|
!     | ~context~     | ~qmckl_context~                        | in     | Global state                |
!     | ~elec_num~    | ~int64_t~                              | in     | Number of electrons         |
!     | ~walk_num~    | ~int64_t~                              | in     | Number of walkers           |
!     | ~coord~       | ~double[3][walk_num][elec_num]~        | in     | Electron coordinates        |
!     | ~ee_distance~ | ~double[walk_num][elec_num][elec_num]~ | out    | Electron-electron distances |


integer function qmckl_compute_ee_distance_f(context, elec_num, walk_num, coord, ee_distance) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: coord(elec_num,walk_num,3)
  double precision      , intent(out) :: ee_distance(elec_num,elec_num,walk_num)

  integer*8 :: k, i, j
  double precision :: x, y, z

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  do k=1,walk_num
     info = qmckl_distance(context, 'T', 'T', elec_num, elec_num, &
          coord(1,k,1), elec_num * walk_num, &
          coord(1,k,1), elec_num * walk_num, &
          ee_distance(1,1,k), elec_num)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
  end do

end function qmckl_compute_ee_distance_f



! #+CALL: generate_c_interface(table=qmckl_ee_distance_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_ee_distance &
    (context, elec_num, walk_num, coord, ee_distance) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: coord(elec_num,3,walk_num)
  real    (c_double ) , intent(out)         :: ee_distance(elec_num,elec_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_ee_distance_f
  info = qmckl_compute_ee_distance_f &
         (context, elec_num, walk_num, coord, ee_distance)

end function qmckl_compute_ee_distance

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_ee_distance_rescaled
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_ee_distance_rescaled_args
!     | Variable                  | Type                                   | In/Out | Description                          |
!     |---------------------------+----------------------------------------+--------+--------------------------------------|
!     | ~context~                 | ~qmckl_context~                        | in     | Global state                         |
!     | ~elec_num~                | ~int64_t~                              | in     | Number of electrons                  |
!     | ~rescale_factor_kappa_ee~ | ~double~                               | in     | Factor to rescale ee distances       |
!     | ~walk_num~                | ~int64_t~                              | in     | Number of walkers                    |
!     | ~coord~                   | ~double[3][walk_num][elec_num]~        | in     | Electron coordinates                 |
!     | ~ee_distance~             | ~double[walk_num][elec_num][elec_num]~ | out    | Electron-electron rescaled distances |


integer function qmckl_compute_ee_distance_rescaled_f(context, elec_num, rescale_factor_kappa_ee, walk_num, &
     coord, ee_distance_rescaled) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: elec_num
  double precision      , intent(in)  :: rescale_factor_kappa_ee
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: coord(elec_num,walk_num,3)
  double precision      , intent(out) :: ee_distance_rescaled(elec_num,elec_num,walk_num)

  integer*8 :: k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  do k=1,walk_num
     info = qmckl_distance_rescaled(context, 'T', 'T', elec_num, elec_num, &
          coord(1,k,1), elec_num * walk_num, &
          coord(1,k,1), elec_num * walk_num, &
          ee_distance_rescaled(1,1,k), elec_num, rescale_factor_kappa_ee)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
  end do

end function qmckl_compute_ee_distance_rescaled_f



! #+CALL: generate_c_interface(table=qmckl_ee_distance_rescaled_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_ee_distance_rescaled &
    (context, elec_num, rescale_factor_kappa_ee, walk_num, coord, ee_distance_rescaled) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_kappa_ee
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: coord(elec_num,3,walk_num)
  real    (c_double ) , intent(out)         :: ee_distance_rescaled(elec_num,elec_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_ee_distance_rescaled_f
  info = qmckl_compute_ee_distance_rescaled_f &
         (context, elec_num, rescale_factor_kappa_ee, walk_num, coord, ee_distance_rescaled)

end function qmckl_compute_ee_distance_rescaled

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_ee_distance_rescaled_deriv_e
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_ee_distance_rescaled_deriv_e_args
!     | Variable                  | Type                                      | In/Out | Description                                     |
!     |---------------------------+-------------------------------------------+--------+-------------------------------------------------|
!     | ~context~                 | ~qmckl_context~                           | in     | Global state                                    |
!     | ~elec_num~                | ~int64_t~                                 | in     | Number of electrons                             |
!     | ~rescale_factor_kappa_ee~ | ~double~                                  | in     | Factor to rescale ee distances                  |
!     | ~walk_num~                | ~int64_t~                                 | in     | Number of walkers                               |
!     | ~coord~                   | ~double[3][walk_num][elec_num]~           | in     | Electron coordinates                            |
!     | ~ee_distance_deriv_e~     | ~double[walk_num][4][elec_num][elec_num]~ | out    | Electron-electron rescaled distance derivatives |


integer function qmckl_compute_ee_distance_rescaled_deriv_e_f(context, elec_num, rescale_factor_kappa_ee, walk_num, &
     coord, ee_distance_rescaled_deriv_e) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: elec_num
  double precision      , intent(in)  :: rescale_factor_kappa_ee
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: coord(elec_num,walk_num,3)
  double precision      , intent(out) :: ee_distance_rescaled_deriv_e(4,elec_num,elec_num,walk_num)

  integer*8 :: k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  do k=1,walk_num
     info = qmckl_distance_rescaled_deriv_e(context, 'T', 'T', elec_num, elec_num, &
          coord(1,k,1), elec_num*walk_num, &
          coord(1,k,1), elec_num*walk_num, &
          ee_distance_rescaled_deriv_e(1,1,1,k), elec_num, rescale_factor_kappa_ee)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
  end do

end function qmckl_compute_ee_distance_rescaled_deriv_e_f



! #+CALL: generate_c_interface(table=qmckl_ee_distance_rescaled_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_ee_distance_rescaled_deriv_e &
    (context, elec_num, rescale_factor_kappa_ee, walk_num, coord, ee_distance_rescaled_deriv_e) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_kappa_ee
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: coord(elec_num,3,walk_num)
  real    (c_double ) , intent(out)         :: ee_distance_rescaled_deriv_e(4,elec_num,elec_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_ee_distance_rescaled_deriv_e_f
  info = qmckl_compute_ee_distance_rescaled_deriv_e_f &
         (context, elec_num, rescale_factor_kappa_ee, walk_num, coord, ee_distance_rescaled_deriv_e)

end function qmckl_compute_ee_distance_rescaled_deriv_e

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_ee_potential
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_ee_potential_args
!     | Variable       | Type                                   | In/Out | Description                          |
!     |----------------+----------------------------------------+--------+--------------------------------------|
!     | ~context~      | ~qmckl_context~                        | in     | Global state                         |
!     | ~elec_num~     | ~int64_t~                              | in     | Number of electrons                  |
!     | ~walk_num~     | ~int64_t~                              | in     | Number of walkers                    |
!     | ~ee_distance~  | ~double[walk_num][elec_num][elec_num]~ | in     | Electron-electron rescaled distances |
!     | ~ee_potential~ | ~double[walk_num]~                     | out    | Electron-electron potential          |


integer function qmckl_compute_ee_potential_f(context, elec_num, walk_num, &
     ee_distance, ee_potential) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: ee_distance(elec_num,elec_num,walk_num)
  double precision      , intent(out) :: ee_potential(walk_num)

  integer*8 :: nw, i, j

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  ee_potential = 0.0d0
  do nw=1,walk_num
    do j=2,elec_num
      do i=1,j-1
        if (dabs(ee_distance(i,j,nw)) > 1e-5) then
          ee_potential(nw) = ee_potential(nw) + 1.0d0/(ee_distance(i,j,nw))
        endif
      end do
    end do
  end do

end function qmckl_compute_ee_potential_f



! #+CALL: generate_c_interface(table=qmckl_ee_potential_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_ee_potential &
    (context, elec_num, walk_num, ee_distance, ee_potential) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: ee_distance(elec_num,elec_num,walk_num)
  real    (c_double ) , intent(out)         :: ee_potential(walk_num)

  integer(c_int32_t), external :: qmckl_compute_ee_potential_f
  info = qmckl_compute_ee_potential_f &
	 (context, elec_num, walk_num, ee_distance, ee_potential)

end function qmckl_compute_ee_potential

! CPU

!    :PROPERTIES:
!    :Name:     qmckl_compute_en_distance
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_en_distance_args
!     | Variable      | Type                                   | In/Out | Description                |
!     |---------------+----------------------------------------+--------+----------------------------|
!     | ~context~     | ~qmckl_context~                        | in     | Global state               |
!     | ~elec_num~    | ~int64_t~                              | in     | Number of electrons        |
!     | ~nucl_num~    | ~int64_t~                              | in     | Number of nuclei           |
!     | ~walk_num~    | ~int64_t~                              | in     | Number of walkers          |
!     | ~elec_coord~  | ~double[3][walk_num][elec_num]~        | in     | Electron coordinates       |
!     | ~nucl_coord~  | ~double[3][elec_num]~                  | in     | Nuclear coordinates        |
!     | ~en_distance~ | ~double[walk_num][nucl_num][elec_num]~ | out    | Electron-nucleus distances |


integer function qmckl_compute_en_distance_f(context, elec_num, nucl_num, walk_num, elec_coord, nucl_coord, en_distance) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: nucl_num
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: elec_coord(elec_num,walk_num,3)
  double precision      , intent(in)  :: nucl_coord(nucl_num,3)
  double precision      , intent(out) :: en_distance(elec_num,nucl_num,walk_num)

  integer*8 :: k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  do k=1,walk_num
     info = qmckl_distance(context, 'T', 'T', elec_num, nucl_num, &
          elec_coord(1,k,1), elec_num * walk_num, &
          nucl_coord, nucl_num, &
          en_distance(1,1,k), elec_num)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
  end do

end function qmckl_compute_en_distance_f



! #+CALL: generate_c_interface(table=qmckl_en_distance_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_en_distance &
    (context, elec_num, nucl_num, walk_num, elec_coord, nucl_coord, en_distance) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: elec_coord(elec_num,walk_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(elec_num,3)
  real    (c_double ) , intent(out)         :: en_distance(elec_num,nucl_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_en_distance_f
  info = qmckl_compute_en_distance_f &
         (context, elec_num, nucl_num, walk_num, elec_coord, nucl_coord, en_distance)

end function qmckl_compute_en_distance

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_en_distance_rescaled
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+NAME: qmckl_en_distance_rescaled_args
!     | Variable                  | Type                                   | In/Out | Description                       |
!     |---------------------------+----------------------------------------+--------+-----------------------------------|
!     | ~context~                 | ~qmckl_context~                        | in     | Global state                      |
!     | ~elec_num~                | ~int64_t~                              | in     | Number of electrons               |
!     | ~nucl_num~                | ~int64_t~                              | in     | Number of nuclei                  |
!     | ~rescale_factor_kappa_en~ | ~double~                               | in     | The factor for rescaled distances |
!     | ~walk_num~                | ~int64_t~                              | in     | Number of walkers                 |
!     | ~elec_coord~              | ~double[3][walk_num][elec_num]~        | in     | Electron coordinates              |
!     | ~nucl_coord~              | ~double[3][elec_num]~                  | in     | Nuclear coordinates               |
!     | ~en_distance_rescaled~    | ~double[walk_num][nucl_num][elec_num]~ | out    | Electron-nucleus distances        |


integer function qmckl_compute_en_distance_rescaled_f(context, elec_num, nucl_num, rescale_factor_kappa_en, walk_num, elec_coord, &
     nucl_coord, en_distance_rescaled) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: nucl_num
  double precision      , intent(in)  :: rescale_factor_kappa_en
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: elec_coord(elec_num,walk_num,3)
  double precision      , intent(in)  :: nucl_coord(nucl_num,3)
  double precision      , intent(out) :: en_distance_rescaled(elec_num,nucl_num,walk_num)

  integer*8 :: k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  ! TODO: comparison with 0
  !if (rescale_factor_kappa_en <= 0) then
  !   info = QMCKL_INVALID_ARG_4
  !   return
  !endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  do k=1,walk_num
     info = qmckl_distance_rescaled(context, 'T', 'T', elec_num, nucl_num, &
          elec_coord(1,k,1), elec_num*walk_num, &
          nucl_coord, nucl_num, &
          en_distance_rescaled(1,1,k), elec_num, rescale_factor_kappa_en)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
  end do

end function qmckl_compute_en_distance_rescaled_f



! #+CALL: generate_c_interface(table=qmckl_en_distance_rescaled_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_en_distance_rescaled &
    (context, elec_num, nucl_num, rescale_factor_kappa_en, walk_num, elec_coord, nucl_coord, en_distance_rescaled) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_kappa_en
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: elec_coord(elec_num,walk_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(elec_num,3)
  real    (c_double ) , intent(out)         :: en_distance_rescaled(elec_num,nucl_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_en_distance_rescaled_f
  info = qmckl_compute_en_distance_rescaled_f &
         (context, elec_num, nucl_num, rescale_factor_kappa_en, walk_num, elec_coord, nucl_coord, en_distance_rescaled)

end function qmckl_compute_en_distance_rescaled

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_en_distance_rescaled_deriv_e
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_en_distance_rescaled_deriv_e_args
!     | Variable                       | Type                                      | In/Out | Description                           |
!     |--------------------------------+-------------------------------------------+--------+---------------------------------------|
!     | ~context~                      | ~qmckl_context~                           | in     | Global state                          |
!     | ~elec_num~                     | ~int64_t~                                 | in     | Number of electrons                   |
!     | ~nucl_num~                     | ~int64_t~                                 | in     | Number of nuclei                      |
!     | ~rescale_factor_kappa_en~      | ~double~                                  | in     | The factor for rescaled distances     |
!     | ~walk_num~                     | ~int64_t~                                 | in     | Number of walkers                     |
!     | ~elec_coord~                   | ~double[3][walk_num][elec_num]~           | in     | Electron coordinates                  |
!     | ~nucl_coord~                   | ~double[3][elec_num]~                     | in     | Nuclear coordinates                   |
!     | ~en_distance_rescaled_deriv_e~ | ~double[walk_num][nucl_num][elec_num][4]~ | out    | Electron-nucleus distance derivatives |


integer function qmckl_compute_en_distance_rescaled_deriv_e_f(context, elec_num, nucl_num, &
     rescale_factor_kappa_en, walk_num, elec_coord, &
     nucl_coord, en_distance_rescaled_deriv_e) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: nucl_num
  double precision      , intent(in)  :: rescale_factor_kappa_en
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: elec_coord(elec_num,walk_num,3)
  double precision      , intent(in)  :: nucl_coord(nucl_num,3)
  double precision      , intent(out) :: en_distance_rescaled_deriv_e(4,elec_num,nucl_num,walk_num)

  integer*8 :: k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  ! TODO: comparison with 0
  !if (rescale_factor_kappa_en <= 0) then
  !   info = QMCKL_INVALID_ARG_4
  !   return
  !endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  do k=1,walk_num
     info = qmckl_distance_rescaled_deriv_e(context, 'T', 'T', elec_num, nucl_num, &
          elec_coord(1,k,1), elec_num*walk_num, &
          nucl_coord, nucl_num, &
          en_distance_rescaled_deriv_e(1,1,1,k), elec_num, rescale_factor_kappa_en)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
  end do

end function qmckl_compute_en_distance_rescaled_deriv_e_f



! #+CALL: generate_c_interface(table=qmckl_en_distance_rescaled_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_en_distance_rescaled_deriv_e &
    (context, elec_num, nucl_num, rescale_factor_kappa_en, walk_num, elec_coord, nucl_coord, en_distance_rescaled_deriv_e) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_kappa_en
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: elec_coord(elec_num,walk_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(elec_num,3)
  real    (c_double ) , intent(out)         :: en_distance_rescaled_deriv_e(elec_num,nucl_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_en_distance_rescaled_deriv_e_f
  info = qmckl_compute_en_distance_rescaled_deriv_e_f &
         (context, elec_num, nucl_num, rescale_factor_kappa_en, walk_num, elec_coord, nucl_coord, en_distance_rescaled_deriv_e)

end function qmckl_compute_en_distance_rescaled_deriv_e

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_en_potential
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_en_potential_args
!     | Variable       | Type                                   | In/Out | Description                          |
!     |----------------+----------------------------------------+--------+--------------------------------------|
!     | ~context~      | ~qmckl_context~                        | in     | Global state                         |
!     | ~elec_num~     | ~int64_t~                              | in     | Number of electrons                  |
!     | ~nucl_num~     | ~int64_t~                              | in     | Number of nuclei                     |
!     | ~walk_num~     | ~int64_t~                              | in     | Number of walkers                    |
!     | ~charge~       | ~double[nucl_num]~                     | in     | charge of nucleus                    |
!     | ~en_distance~  | ~double[walk_num][nucl_num][elec_num]~ | in     | Electron-electron rescaled distances |
!     | ~en_potential~ | ~double[walk_num]~                     | out    | Electron-electron potential          |


integer function qmckl_compute_en_potential_f(context, elec_num, nucl_num, walk_num, &
     charge, en_distance, en_potential) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: elec_num
  integer*8             , intent(in)  :: nucl_num
  integer*8             , intent(in)  :: walk_num
  double precision      , intent(in)  :: charge(nucl_num)
  double precision      , intent(in)  :: en_distance(elec_num,nucl_num,walk_num)
  double precision      , intent(out) :: en_potential(walk_num)

  integer*8 :: nw, i, j

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  en_potential = 0.0d0
  do nw=1,walk_num
    do j=1,nucl_num
      do i=1,elec_num
        if (dabs(en_distance(i,j,nw)) > 1e-5) then
          en_potential(nw) = en_potential(nw) - charge(j)/(en_distance(i,j,nw))
        endif
      end do
    end do
  end do

end function qmckl_compute_en_potential_f



! #+CALL: generate_c_interface(table=qmckl_en_potential_args,rettyp=get_value("CRetType"),fname=get_value("Name"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_en_potential &
    (context, elec_num, nucl_num, walk_num, charge, en_distance, en_potential) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: charge(nucl_num)
  real    (c_double ) , intent(in)          :: en_distance(elec_num,nucl_num,walk_num)
  real    (c_double ) , intent(out)         :: en_potential(walk_num)

  integer(c_int32_t), external :: qmckl_compute_en_potential_f
  info = qmckl_compute_en_potential_f &
	 (context, elec_num, nucl_num, walk_num, charge, en_distance, en_potential)

end function qmckl_compute_en_potential
