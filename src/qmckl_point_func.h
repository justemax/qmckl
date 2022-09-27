/* Number of points */


qmckl_exit_code qmckl_get_point_num (const qmckl_context context, int64_t* const num);

/* Point coordinates */


qmckl_exit_code qmckl_get_point(const qmckl_context context,
                                const char transp,
                                double* const coord,
                                const int64_t size_max);

/* Initialization functions */

/*    When the data is set in the context, if the arrays are large */
/*    enough, we overwrite the data contained in them. */

/*    To set the data relative to the points in the context, one of the */
/*    following functions need to be called. Here, ~num~ is the number of */
/*    points to set. */


qmckl_exit_code qmckl_set_point (qmckl_context context,
                                 const char transp,
                                 const int64_t num,
                                 const double* coord,
                                 const int64_t size_max);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_set_point_device (qmckl_context context,
                                        const char transp,
                                        const int64_t num,
                                        const double* coord,
										const int64_t size_max,
                                        int device_id);
#endif
