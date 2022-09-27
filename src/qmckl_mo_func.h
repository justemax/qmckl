/* Initialization functions */

/*    To set the basis set, all the following functions need to be */
/*    called. */


qmckl_exit_code  qmckl_set_mo_basis_mo_num           (qmckl_context context, const int64_t   mo_num);
qmckl_exit_code  qmckl_set_mo_basis_coefficient      (qmckl_context context, const double  * coefficient);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code  qmckl_set_mo_basis_coefficient_device (qmckl_context context, const double  * coefficient, int device_id);
#endif

/* Access functions */


qmckl_exit_code
qmckl_get_mo_basis_mo_num (const qmckl_context context,
                           int64_t* mo_num);

qmckl_exit_code
qmckl_get_mo_basis_coefficient (const qmckl_context context,
                                double* const coefficient,
                                const int64_t size_max);



/* When all the data for the AOs have been provided, the following */
/* function returns ~true~. */


bool qmckl_mo_basis_provided (const qmckl_context context);

/* Update */

/*    Useless MOs can be removed, for instance virtual MOs in a single */
/*    determinant calculation. */

/*    To select a subset of MOs that will be kept, create an array of */
/*    integers of size =mo_num=. If the integer is zero, the MO is dropped, */
/*    otherwise it is kept. */


bool qmckl_mo_basis_select_mo (const qmckl_context context,
                               const int32_t* keep,
                               const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_mo_basis_mo_value(qmckl_context context,
                            double* const mo_value,
                            const int64_t size_max);



/* Uses the given array to compute the values. */


qmckl_exit_code
qmckl_get_mo_basis_mo_value_inplace (qmckl_context context,
                                     double* const mo_value,
                                     const int64_t size_max);



/*  #+CALL: generate_c_header(table=qmckl_mo_basis_mo_value_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_value (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const double* coefficient_t,
      const double* ao_value,
      double* const mo_value );



/* #+CALL: generate_c_header(table=qmckl_mo_basis_mo_value_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value_doc")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_value_doc (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const double* coefficient_t,
      const double* ao_value,
      double* const mo_value );

/* HPC version */



#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_value_hpc (const qmckl_context context,
                                     const int64_t ao_num,
                                     const int64_t mo_num,
                                     const int64_t point_num,
                                     const double* coefficient_t,
                                     const double* ao_value,
                                     double* const mo_value );
#endif

/* Get */


qmckl_exit_code
qmckl_get_mo_basis_mo_vgl(qmckl_context context,
                          double* const mo_vgl,
                          const int64_t size_max);



/* Uses the given array to compute the VGL. */


qmckl_exit_code
qmckl_get_mo_basis_mo_vgl_inplace (qmckl_context context,
                                   double* const mo_vgl,
                                   const int64_t size_max);



/*  #+CALL: generate_c_header(table=qmckl_mo_basis_mo_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const double* coefficient_t,
      const double* ao_vgl,
      double* const mo_vgl );



/* #+CALL: generate_c_header(table=qmckl_mo_basis_mo_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl_doc")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_doc (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const double* coefficient_t,
      const double* ao_vgl,
      double* const mo_vgl );

/* HPC version */



#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl_hpc (const qmckl_context context,
                                   const int64_t ao_num,
                                   const int64_t mo_num,
                                   const int64_t point_num,
                                   const double* coefficient_t,
                                   const double* ao_vgl,
                                   double* const mo_vgl );
#endif
