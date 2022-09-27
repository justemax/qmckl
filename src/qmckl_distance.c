#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include <stdio.h>

#include "qmckl.h"

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_distance_device (
    const qmckl_context context,
    const char transa,
    const char transb,
    const int64_t m,
    const int64_t n,
    const double* A,
    const int64_t lda,
    const double* B,
    const int64_t ldb,
    double* const C,
    const int64_t ldc,
    int device_id
) {
  uint64_t transab;
  double x,y,z;

  if (context == QMCKL_NULL_CONTEXT){
    return QMCKL_INVALID_CONTEXT;
  }

  if(m <= 1e-5){
    return QMCKL_INVALID_ARG_4;
  }

  if(n <= 1e-5){
    return QMCKL_INVALID_ARG_5;
  }


  if (transa == 'N' || transa == 'n'){
    transab = 0;
  }
  else if (transa == 'T' || transa == 't'){
    transab = 1;
  }else{
    transab = -100;
  }

  if (transb == 'N' || transb == 'n'){
    transab = transab;
  }
  else if (transa == 'T' || transa == 't'){
    transab = transab + 2;
  }else{
    transab = -100;
  }

  // check for LDA
  if (transab & 1 == 0 && lda < 3) {
    return QMCKL_INVALID_ARG_7;
  }

  if (transab & 1 == 1 && lda < m) {
    return QMCKL_INVALID_ARG_7;
  }

  if (transab & 2 == 0 && lda < 3) {
    return QMCKL_INVALID_ARG_7;
  }

  if (transab & 2 == 2 && lda < m) {
    return QMCKL_INVALID_ARG_7;
  }

  // check for LDB
  if (transab & 1 == 0 && ldb < 3) {
    return QMCKL_INVALID_ARG_9;
  }

  if (transab & 1 == 1 && ldb < n) {
    return QMCKL_INVALID_ARG_9;
  }

  if (transab & 2 == 0 && ldb < 3) {
    return QMCKL_INVALID_ARG_9;
  }

  if (transab & 2 == 2 && ldb < n) {
    return QMCKL_INVALID_ARG_9;
  }
  // check for LDC
  if (ldc < m) {
    return QMCKL_INVALID_ARG_11;
  }

  switch(transab)
    {
    case(0):
      #pragma omp target is_device_ptr(A, B, C)
      {
      #pragma omp teams distribute parallel for simd
      for(int j = 0; j < n; j++){
        for(int i = 0; j < m; i++){
          x = A[i * 3 + 0] - B[j * 3 + 0];
          y = A[i * 3 + 1] - B[j * 3 + 1];
          z = A[i * 3 + 2] - B[j * 3 + 2];

          C[j * n + i] = sqrt(x*x + y*y + z*z);
        }
      }
      }
      break;

    case(1):
     #pragma omp target is_device_ptr(A, B, C)
     {
     #pragma omp teams distribute parallel for simd
     for(int j = 0; j < n; j++){
        for(int i = 0; i < m; i++){
          x = A[0 * lda + i] - B[j * 3 + 0];
          y = A[1 * lda + i] - B[j * 3 + 1];
          z = A[2 * lda + i] - B[j * 3 + 2];

          C[j * n + i] = sqrt(x*x + y*y + z*z);
        }
      }
      }
      break;

    case(2):
      #pragma omp target is_device_ptr(A, B, C)
      {
      #pragma omp teams distribute parallel for simd
      for(int j = 0; j < n; j++){
        for(int i = 0; i < m; i++){
          x = A[i * 3 + 0] - B[0 * ldb + j];
          y = A[i * 3 + 1] - B[1 * ldb + j];
          z = A[i * 3 + 2] - B[2 * ldb + j];

          C[j * n + i] = sqrt(x*x + y*y + z*z);
        }
      }
      }
      break;

    case(3):
      #pragma omp target is_device_ptr(A, B, C)
      {
      #pragma omp teams distribute parallel for simd
      for(int j = 0; j < n; j++){
        for(int i = 0; i < m; i++){
          x = A[0 * lda + i] - B[0 * ldb + j];
          y = A[1 * lda + i] - B[1 * ldb + j];
          z = A[2 * lda + i] - B[2 * ldb + j];

          C[j * n + i] = sqrt((x*x) + (y*y) + (z*z));
        }
      }
      }
      break;
    }


  return QMCKL_SUCCESS;

}
#endif
