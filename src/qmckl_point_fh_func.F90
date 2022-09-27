interface
  integer(c_int32_t) function qmckl_get_point_num(context, num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_point(context, transp, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double ) , intent(out)         :: coord(*)
    integer (c_int64_t) , intent(in)          :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_point(context, &
       transp, num, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    integer (c_int64_t) , intent(in)  , value :: num
    real    (c_double ) , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

#ifdef HAVE_DEVICE_POINTERS
interface
  integer(c_int32_t) function qmckl_set_point_device(context, &
       transp, num, coord, size_max, device_id) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    integer (c_int64_t) , intent(in)  , value :: num
    real    (c_double ) , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
    integer (c_int32_t) , intent(in)  , value :: device_id
end interface
#endif
