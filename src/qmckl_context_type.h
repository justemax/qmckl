/* Context handling */

/*   The context variable is a handle for the state of the library, */
/*   and is stored in a data structure which can't be seen outside of */
/*   the library. To simplify compatibility with other languages, the */
/*   pointer to the internal data structure is converted into a 64-bit */
/*   signed integer, defined in the ~qmckl_context~ type. */
/*   A value of ~QMCKL_NULL_CONTEXT~ for the context is equivalent to a */
/*   ~NULL~ pointer. */

/*   #+NAME: qmckl_context */

typedef int64_t qmckl_context ;
#define QMCKL_NULL_CONTEXT (qmckl_context) 0
