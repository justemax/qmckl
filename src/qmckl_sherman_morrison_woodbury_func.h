/* C header */

/*     #+CALL: generate_c_header(table=qmckl_sherman_morrison_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_sherman_morrison(
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant);

/* C header */

/*     #+CALL: generate_c_header(table=qmckl_woodbury_2_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_woodbury_2(
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant);

/* C header */

/*     #+CALL: generate_c_header(table=qmckl_woodbury_3_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_woodbury_3(
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
double* determinant);

/* C header */

/*     #+CALL: generate_c_header(table=qmckl_sherman_morrison_splitting_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_sherman_morrison_splitting(
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant);

/* C header */

/*     #+CALL: generate_c_header(table=qmckl_sherman_morrison_smw32s_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_sherman_morrison_smw32s(
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant);

/* C header */

/*     #+CALL: generate_c_header(table=qmckl_slagel_splitting_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_slagel_splitting (
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* later_updates,
      uint64_t* later_index,
      uint64_t* later,
      double* determinant);
