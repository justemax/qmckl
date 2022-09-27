


/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */


qmckl_exit_code qmckl_init_jastrow(qmckl_context context);

/* Access functions */


qmckl_exit_code  qmckl_get_jastrow_aord_num          (qmckl_context context, int64_t* const aord_num);
qmckl_exit_code  qmckl_get_jastrow_bord_num          (qmckl_context context, int64_t* const bord_num);
qmckl_exit_code  qmckl_get_jastrow_cord_num          (qmckl_context context, int64_t* const bord_num);
qmckl_exit_code  qmckl_get_jastrow_type_nucl_num     (qmckl_context context, int64_t* const type_nucl_num);
qmckl_exit_code  qmckl_get_jastrow_type_nucl_vector  (qmckl_context context, int64_t* const type_nucl_num, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_aord_vector       (qmckl_context context, double * const aord_vector, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_bord_vector       (qmckl_context context, double * const bord_vector, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_cord_vector       (qmckl_context context, double * const cord_vector, const int64_t size_max);



/* Along with these core functions, calculation of the jastrow factor */
/* requires the following additional information to be set: */


/* When all the data for the AOs have been provided, the following */
/* function returns ~true~. */


bool      qmckl_jastrow_provided           (const qmckl_context context);

/* Host initialization */

/*    To prepare for the Jastrow and its derivative, all the following functions need to be */
/*    called. */


qmckl_exit_code  qmckl_set_jastrow_ord_num           (qmckl_context context, const int64_t aord_num, const int64_t bord_num, const int64_t cord_num);
qmckl_exit_code  qmckl_set_jastrow_type_nucl_num     (qmckl_context context, const int64_t type_nucl_num);
qmckl_exit_code  qmckl_set_jastrow_type_nucl_vector  (qmckl_context context, const int64_t* type_nucl_vector, const int64_t nucl_num);
qmckl_exit_code  qmckl_set_jastrow_aord_vector       (qmckl_context context, const double * aord_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_bord_vector       (qmckl_context context, const double * bord_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_cord_vector       (qmckl_context context, const double * cord_vector, const int64_t size_max);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code  qmckl_set_jastrow_type_nucl_vector_device  (qmckl_context context, const int64_t* type_nucl_vector, const int64_t nucl_num, int device_id);
qmckl_exit_code  qmckl_set_jastrow_aord_vector_device       (qmckl_context context, const double * aord_vector, const int64_t size_max, int device_id);
qmckl_exit_code  qmckl_set_jastrow_bord_vector_device       (qmckl_context context, const double * bord_vector, const int64_t size_max, int device_id);
qmckl_exit_code  qmckl_set_jastrow_cord_vector_device       (qmckl_context context, const double * cord_vector, const int64_t size_max, int device_id);
#endif

/* Get */

qmckl_exit_code
qmckl_get_jastrow_asymp_jasb(qmckl_context context,
                             double* const asymp_jasb,
                             const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_factor_ee(qmckl_context context,
                            double* const factor_ee,
                            const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_factor_ee_deriv_e(qmckl_context context,
                                    double* const factor_ee_deriv_e,
                                    const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_factor_en(qmckl_context context,
                            double* const factor_en,
                            const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_factor_en_deriv_e(qmckl_context context,
                                    double* const factor_en_deriv_e,
                                    const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_een_rescaled_e(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_een_rescaled_e_deriv_e(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_een_rescaled_n(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_een_rescaled_n_deriv_e(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max);

/* CPU & offload version */


qmckl_exit_code qmckl_get_jastrow_dim_cord_vect(qmckl_context context, int64_t* const dim_cord_vect);
qmckl_exit_code qmckl_get_jastrow_cord_vect_full(qmckl_context context, double* const cord_vect_full);
qmckl_exit_code qmckl_get_jastrow_lkpm_combined_index(qmckl_context context, int64_t* const lkpm_combined_index);
qmckl_exit_code qmckl_get_jastrow_tmp_c(qmckl_context context, double* const tmp_c);
qmckl_exit_code qmckl_get_jastrow_dtmp_c(qmckl_context context, double* const dtmp_c);

/* Device pointers version */


qmckl_exit_code qmckl_get_jastrow_cord_vect_full_device(qmckl_context context, double* const cord_vect_full, int device_id);
qmckl_exit_code qmckl_get_jastrow_lkpm_combined_index_device(qmckl_context context, int64_t* const lkpm_combined_index, int device_id);
qmckl_exit_code qmckl_get_jastrow_tmp_c_device(qmckl_context context, double* const tmp_c, int device_id);
qmckl_exit_code qmckl_get_jastrow_dtmp_c_device(qmckl_context context, double* const dtmp_c, int device_id);





/* #+CALL: generate_c_header(table=qmckl_factor_tmp_c_args,rettyp=get_value("CRetType"),fname="qmckl_compute_tmp_c") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_tmp_c (
          const qmckl_context context,
          const int64_t cord_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* een_rescaled_e,
          const double* een_rescaled_n,
          double* const tmp_c );

/* Get */

qmckl_exit_code
qmckl_get_jastrow_factor_een(qmckl_context context,
                             double* const factor_een,
                             const int64_t size_max);

/* CPU */


qmckl_exit_code
qmckl_get_jastrow_factor_een_deriv_e(qmckl_context context,
                                     double* const factor_een_deriv_e,
                                     const int64_t size_max);

/* Device pointers */


#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_get_jastrow_factor_een_deriv_e_device(
  qmckl_context context,
  double* const factor_een_deriv_e,
  const int64_t size_max,
  int device_id
);
#endif
