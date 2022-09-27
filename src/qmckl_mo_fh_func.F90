! Fortran interfaces


interface
  integer(c_int32_t) function qmckl_get_mo_basis_mo_num (context, &
       mo_num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: mo_num
  end function qmckl_get_mo_basis_mo_num
end interface

interface
  integer(c_int32_t) function qmckl_get_mo_basis_coefficient(context, &
       coefficient, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    double precision, intent(out)             :: coefficient(*)
    integer (c_int64_t) , intent(in), value   :: size_max
  end function qmckl_get_mo_basis_coefficient
end interface

! Fortran interface


interface
  integer(c_int32_t) function qmckl_mo_basis_select_mo (context, &
       keep, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in), value :: context
    integer (c_int32_t) , intent(in)        :: keep(*)
    integer (c_int64_t) , intent(in), value :: size_max
  end function qmckl_mo_basis_select_mo
end interface

interface
   integer(c_int32_t) function qmckl_get_mo_basis_mo_value (context, &
        mo_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: mo_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_value
end interface

interface
   integer(c_int32_t) function qmckl_get_mo_basis_mo_value_inplace (context, &
        mo_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: mo_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_value_inplace
end interface

interface
   integer(c_int32_t) function qmckl_get_mo_basis_mo_vgl (context, &
        mo_vgl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: mo_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_vgl
end interface

interface
   integer(c_int32_t) function qmckl_get_mo_basis_mo_vgl_inplace (context, &
        mo_vgl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: mo_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_vgl_inplace
end interface
