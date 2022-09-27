interface
   subroutine qmckl_string_of_error (error, string) bind(C, name='qmckl_string_of_error_f')
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_exit_code), intent(in), value :: error
     character, intent(out) :: string(128)
   end subroutine qmckl_string_of_error
end interface

interface
   subroutine qmckl_last_error (context, string) bind(C, name='qmckl_last_error')
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in), value :: context
     character, intent(out) :: string(*)
   end subroutine qmckl_last_error
end interface
