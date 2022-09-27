/* C source */


#include <stdbool.h>
#include <math.h>
#include "qmckl.h"

qmckl_exit_code qmckl_sherman_morrison(const qmckl_context context,
                                const uint64_t LDS,
                                const uint64_t Dim,
                                const uint64_t N_updates,
                                const double* Updates,
                                const uint64_t* Updates_index,
                                const double breakdown,
                                double* Slater_inv,
                                double* determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  double C[Dim];
  double D[Dim];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < Dim; i++) {
      C[i] = 0;
      for (uint64_t j = 0; j < Dim; j++) {
        C[i] += Slater_inv[i * LDS + j] * Updates[l * Dim + j];
      }
    }

    // Denominator
    double den = 1 + C[Updates_index[l] - 1];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1 / den;

    // Update det(A)
    if (determinant != NULL)
      *determinant *= den;

    // D = v^T x A^{-1}
    for (uint64_t j = 0; j < Dim; j++) {
      D[j] = Slater_inv[(Updates_index[l] - 1) * LDS + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < Dim; i++) {
      for (uint64_t j = 0; j < Dim; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * LDS + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

/* C source */


#include <stdbool.h>
#include <math.h>
#include "qmckl.h"

qmckl_exit_code qmckl_woodbury_2(const qmckl_context context,
                                const uint64_t LDS,
                                const uint64_t Dim,
                                const double* Updates,
                                const uint64_t* Updates_index,
                                const double breakdown,
                                double* Slater_inv,
                                double* determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = S_inv * U  !! NON-STANDARD MATRIX MULTIPLICATION BECAUSE
  // OF LAYOUT OF 'Updates' !!
  double C[2 * Dim];
  for (uint64_t i = 0; i < Dim; i++) {
    for (uint64_t j = 0; j < 2; j++) {
      C[i * 2 + j] = 0;
      for (uint64_t k = 0; k < Dim; k++) {
        C[i * 2 + j] += Slater_inv[i * LDS + k] * Updates[Dim * j + k];
      }
    }
  }

  // Compute B = 1 + V * C
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
    return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant != NULL)
    *determinant *= det;

  // Compute B^{-1} with explicit formula for 2x2 inversion
  double Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // Compute tmp = B^{-1} x (V.S^{-1})
  double tmp[2 * Dim];
  for (uint64_t i = 0; i < 2; i++) {
    for (uint64_t j = 0; j < Dim; j++) {
      tmp[i * Dim + j] = Binv[i * 2] * Slater_inv[row1 * LDS + j];
      tmp[i * Dim + j] += Binv[i * 2 + 1] * Slater_inv[row2 * LDS + j];
    }
  }

  // Compute (S + U V)^{-1} = S^{-1} - C x tmp
  for (uint64_t i = 0; i < Dim; i++) {
    for (uint64_t j = 0; j < Dim; j++) {
      Slater_inv[i * LDS + j] -= C[i * 2] * tmp[j];
      Slater_inv[i * LDS + j] -= C[i * 2 + 1] * tmp[Dim + j];
    }
  }

  return QMCKL_SUCCESS;
}

/* C source */


#include <stdbool.h>
#include <math.h>
#include "qmckl.h"

qmckl_exit_code qmckl_woodbury_3(const qmckl_context context,
                                const uint64_t LDS,
                                const uint64_t Dim,
                                const double* Updates,
                                const uint64_t* Updates_index,
                                const double breakdown,
                                double* Slater_inv,
                                double* determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = S_inv * U  !! NON-STANDARD MATRIX MULTIPLICATION BECAUSE
  // OF LAYOUT OF 'Updates' !!
  double C[3 * Dim];
  for (uint64_t i = 0; i < Dim; i++) {
    for (uint64_t j = 0; j < 3; j++) {
      C[i * 3 + j] = 0;
      for (uint64_t k = 0; k < Dim; k++) {
        C[i * 3 + j] += Slater_inv[i * LDS + k] * Updates[Dim * j + k];
      }
    }
  }

  // Compute B = 1 + V.C
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
    return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant != NULL)
    *determinant *= det;

  // Compute B^{-1} with explicit formula for 3x3 inversion
  double Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // Compute tmp = B^{-1} x (V.S^{-1})
  double tmp[3 * Dim];
  for (uint64_t i = 0; i < 3; i++) {
    for (uint64_t j = 0; j < Dim; j++) {
      tmp[i * Dim + j] = Binv[i * 3] * Slater_inv[row1 * LDS + j];
      tmp[i * Dim + j] += Binv[i * 3 + 1] * Slater_inv[row2 * LDS + j];
      tmp[i * Dim + j] += Binv[i * 3 + 2] * Slater_inv[row3 * LDS + j];
    }
  }

  // Compute (S + U V)^{-1} = S^{-1} - C x tmp
  for (uint64_t i = 0; i < Dim; i++) {
    for (uint64_t j = 0; j < Dim; j++) {
      Slater_inv[i * LDS + j] -= C[i * 3] * tmp[j];
      Slater_inv[i * LDS + j] -= C[i * 3 + 1] * tmp[Dim + j];
      Slater_inv[i * LDS + j] -= C[i * 3 + 2] * tmp[2 * Dim + j];
    }
  }

  return QMCKL_SUCCESS;
}

/* C source */


#include <stdbool.h>
#include "qmckl.h"

qmckl_exit_code qmckl_sherman_morrison_splitting(const qmckl_context context,
                                const uint64_t LDS,
                                const uint64_t Dim,
                                const uint64_t N_updates,
                                const double* Updates,
                                const uint64_t* Updates_index,
                                const double breakdown,
                                double* Slater_inv,
                                double* determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  double later_updates[Dim * N_updates];
  uint64_t later_index[N_updates];
  uint64_t later = 0;

  (void) qmckl_slagel_splitting(LDS, Dim, N_updates, Updates, Updates_index,
                            breakdown, Slater_inv, later_updates, later_index, &later, determinant);

  if (later > 0) {
    (void) qmckl_sherman_morrison_splitting(context, LDS, Dim, later,
                            later_updates, later_index, breakdown, Slater_inv, determinant);
  }

  return QMCKL_SUCCESS;
}

/* C source */


#include <stdbool.h>
#include "qmckl.h"

qmckl_exit_code qmckl_sherman_morrison_smw32s(const qmckl_context context,
                                const uint64_t LDS,
                                const uint64_t Dim,
                                const uint64_t N_updates,
                                const double* Updates,
                                const uint64_t* Updates_index,
                                const double breakdown,
                                double* Slater_inv,
                                double* determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  uint64_t n_of_3blocks = N_updates / 3;
  uint64_t remainder = N_updates % 3;
  uint64_t length_3block = 3 * Dim;

  // Apply first 3*n_of_3blocks updates in n_of_3blocks blocks of 3 updates with
  // Woodbury 3x3 kernel
  double later_updates[Dim * N_updates];
  uint64_t later_index[N_updates];
  uint64_t later = 0;
  if (n_of_3blocks > 0) {
    for (uint64_t i = 0; i < n_of_3blocks; i++) {
      const double *Updates_3block = &Updates[i * length_3block];
      const uint64_t *Updates_index_3block = &Updates_index[i * 3];
      rc = qmckl_woodbury_3(context, LDS, Dim, Updates_3block, Updates_index_3block, breakdown, Slater_inv, determinant);
      if (rc != 0) { // Send the entire block to slagel_splitting
        uint64_t l = 0;
        (void) qmckl_slagel_splitting(LDS, Dim, 3, Updates_3block, Updates_index_3block,
                breakdown, Slater_inv, later_updates + (Dim * later), later_index + later, &l, determinant);
        later = later + l;
      }
    }
  }

  // Apply last remaining block of 2 updates with Woodbury 2x2 kernel
  if (remainder == 2) {
    const double *Updates_2block = &Updates[n_of_3blocks * length_3block];
    const uint64_t *Updates_index_2block = &Updates_index[3 * n_of_3blocks];
    rc = qmckl_woodbury_2(context, LDS, Dim, Updates_2block, Updates_index_2block, breakdown, Slater_inv, determinant);
    if (rc != 0) { // Send the entire block to slagel_splitting
      uint64_t l = 0;
      (void) qmckl_slagel_splitting(LDS, Dim, 2, Updates_2block, Updates_index_2block,
              breakdown, Slater_inv, later_updates + (Dim * later), later_index + later, &l, determinant);
      later = later + l;
    }
  }
  // Apply last remaining update with slagel_splitting
  else if (remainder == 1) {
    const double *Updates_1block = &Updates[n_of_3blocks * length_3block];
    const uint64_t *Updates_index_1block = &Updates_index[3 * n_of_3blocks];
    uint64_t l = 0;
    (void) qmckl_slagel_splitting(LDS, Dim, 1, Updates_1block, Updates_index_1block,
            breakdown, Slater_inv, later_updates + (Dim * later), later_index + later, &l, determinant);
    later = later + l;
  }

  if (later > 0) {
    (void) qmckl_sherman_morrison_splitting(context, LDS, Dim, later, later_updates, later_index, breakdown, Slater_inv, determinant);
  }
  return QMCKL_SUCCESS;
}

/* C source */


#include <stdbool.h>
#include <math.h>
#include "qmckl.h"

qmckl_exit_code qmckl_slagel_splitting(uint64_t LDS,
                                              uint64_t Dim,
                                              uint64_t N_updates,
                                              const double *Updates,
                                              const uint64_t *Updates_index,
                                              const double breakdown,
                                              double *Slater_inv,
                                              double *later_updates,
                                              uint64_t *later_index,
                                              uint64_t *later,
                                              double *determinant) {
// #ifdef DEBUG //  Leave commented out since debugging information is not yet implemented in QMCkl.
//   std::cerr << "Called slagel_splitting with " << N_updates << " updates" << std::endl;
// #endif

  double C[Dim];
  double D[Dim];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < Dim; i++) {
      C[i] = 0;
      for (uint64_t j = 0; j < Dim; j++) {
        C[i] += Slater_inv[i * LDS + j] * Updates[l * Dim + j];
      }
    }

    // Denominator
    double den = 1 + C[Updates_index[l] - 1];
    if (fabs(den) < breakdown) { // Here is decided to split the update, or not.

      // U_l = U_l / 2: split the update in 2 equal halves and save the second halve in later_updates
      for (uint64_t i = 0; i < Dim; i++) {
        later_updates[*later * Dim + i] = Updates[l * Dim + i] / 2.0;
        C[i] /= 2.0;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1 + C[Updates_index[l] - 1];
    } // From here onwards we continue with applying the first havel of the update to Slater_inv
    double iden = 1 / den;

    if (determinant != NULL)
      *determinant *= den;

    // D = v^T x S^{-1}
    for (uint64_t j = 0; j < Dim; j++) {
      D[j] = Slater_inv[(Updates_index[l] - 1) * LDS + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < Dim; i++) {
      for (uint64_t j = 0; j < Dim; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * LDS + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}
