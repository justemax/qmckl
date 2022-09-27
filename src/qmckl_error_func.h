/* Decoding errors */

/*    To decode the error messages, ~qmckl_string_of_error~ converts an */
/*    error code into a string. */


const char*
qmckl_string_of_error (const qmckl_exit_code error);

/* Updating errors in the context */

/*    The error is updated in the context using ~qmckl_set_error~. */
/*    When the error is set in the context, it is mandatory to specify */
/*    from which function the error is triggered, and a message */
/*    explaining the error. The exit code can't be ~QMCKL_SUCCESS~. */

/*    # Header */

qmckl_exit_code
qmckl_set_error(qmckl_context context,
                const qmckl_exit_code exit_code,
                const char* function_name,
                const char* message);

/* Get the error */

/*   Upon error, the error type and message can be obtained from the */
/*   context using ~qmckl_get_error~. The message and function name */
/*   is returned in the variables provided. Therefore, passing a  */
/*   function name and message is mandatory. */

/*   # Header */

qmckl_exit_code
qmckl_get_error(qmckl_context context,
                qmckl_exit_code *exit_code,
                char* function_name,
                char* message);

/* Failing */
  
/*    To make a function fail, the ~qmckl_failwith~ function should be */
/*    called, such that information about the failure is stored in */
/*    the context. The desired exit code is given as an argument, as */
/*    well as the name of the function and an error message. If the */
/*    message is ~NULL~, then the default message obtained by */
/*    ~qmckl_string_of_error~ is used. The return code of the function is */
/*    the desired return code. */
/*    Upon failure, a ~QMCKL_NULL_CONTEXT~ is returned. */


qmckl_exit_code
qmckl_failwith(qmckl_context context,
               const qmckl_exit_code exit_code,
               const char* function,
               const char* message) ;

/* Last error */

/*   Returns a string describing the last error, using ~qmckl_get_error~. */

/*   # Header */

qmckl_exit_code
qmckl_last_error(qmckl_context context, char* buffer);
