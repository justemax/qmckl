qmckl_exit_code
qmckl_get_nucleus_num(const qmckl_context context,
                      int64_t* const num);

qmckl_exit_code
qmckl_get_nucleus_charge(const qmckl_context context,
                         double* const charge,
                         const int64_t size_max);

qmckl_exit_code
qmckl_get_nucleus_rescale_factor(const qmckl_context context,
                                 double* const rescale_factor_kappa);

qmckl_exit_code
qmckl_get_nucleus_coord(const qmckl_context context,
                        const char transp,
                        double* const coord,
                        const int64_t size_max);



/* When all the data relative to nuclei have been set, the following */
/* function returns ~true~. */


bool qmckl_nucleus_provided (const qmckl_context context);



/* To set the data relative to the nuclei in the context, the */
/* following functions need to be called. */


qmckl_exit_code
qmckl_set_nucleus_num(qmckl_context context,
                      const int64_t num);

qmckl_exit_code
qmckl_set_nucleus_charge(qmckl_context context,
                         const double* charge,
                         const int64_t size_max);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_set_nucleus_charge_device(qmckl_context context,
                                const double* charge,
                                const int64_t size_max,
								int device_id);
#endif

qmckl_exit_code
qmckl_set_nucleus_coord(qmckl_context context,
                        const char transp,
                        const double* coord,
                        const int64_t size_max);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_set_nucleus_coord_device(qmckl_context context,
                               const char transp,
                               const double* coord,
                               const int64_t size_max,
							   int device_id);
#endif

qmckl_exit_code
qmckl_set_nucleus_rescale_factor(qmckl_context context,
                                 const double kappa);

/* Get */


qmckl_exit_code
qmckl_get_nucleus_nn_distance(qmckl_context context,
                              double* distance,
                              const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_nucleus_nn_distance_rescaled(qmckl_context context,
                                       double* distance_rescaled,
                                       const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_nucleus_repulsion(qmckl_context context, double* const energy);
