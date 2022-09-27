interface
  integer(c_int32_t) function qmckl_get_nucleus_num(context, num) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function qmckl_get_nucleus_num
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_charge(context, charge, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: charge(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_nucleus_charge
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_rescale_factor(context, kappa) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: kappa
  end function qmckl_get_nucleus_rescale_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_coord(context, transp, coord, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double)  , intent(out)         :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_nucleus_coord
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_num(context, num) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_nucleus_num
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_charge(context, charge, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: charge(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_coord(context, transp, coord, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double)  , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_nucleus_coord
end interface

#ifdef HAVE_DEVICE_POINTERS
interface
  integer(c_int32_t) function qmckl_set_nucleus_coord_device(context, transp, coord, size_max, device_id) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double)  , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
	integer (c_int32_t) , intent(in)  , value :: device_id
  end function qmckl_set_nucleus_coord_device
end interface
#endif

interface
  integer(c_int32_t) function qmckl_set_nucleus_rescale_factor(context, kappa) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)  , value :: kappa
  end function qmckl_set_nucleus_rescale_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_nn_distance(context, distance, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_nn_distance_rescaled(context, distance_rescaled, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance_rescaled(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_repulsion(context, energy) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: energy
  end function
end interface
