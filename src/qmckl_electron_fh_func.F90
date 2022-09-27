interface
  integer(c_int32_t) function qmckl_set_electron_num(context, alpha, beta) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: beta
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_electron_coord(context, transp, walk_num, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transp
    integer (c_int64_t) , intent(in)  , value :: walk_num
    double precision    , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_electron_coord_device(context, transp, walk_num, coord, size_max, device_id) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transp
    integer (c_int64_t) , intent(in)  , value :: walk_num
    double precision    , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
    integer (c_int32_t) , intent(in)  , value :: device_id
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_electron_ee_distance(context, distance) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_electron_en_distance(context, distance) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
  end function
end interface
