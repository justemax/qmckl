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
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_jastrow_private_func.h"
#include "qmckl_jastrow_private_type.h"

#ifdef HAVE_CUBLAS_OFFLOAD
#include "cublas_v2.h"
#endif

#ifdef HAVE_DEVICE_POINTERS
#include <omp.h>
#endif

qmckl_exit_code qmckl_init_jastrow(qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  ctx->jastrow.uninitialized = (1 << 6) - 1;

  /* Default values */
  return QMCKL_SUCCESS;
}

bool qmckl_jastrow_provided(const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  return ctx->jastrow.provided;
}

qmckl_exit_code qmckl_get_jastrow_aord_num (const qmckl_context context, int64_t* const aord_num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  if (aord_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_aord_num",
                           "aord_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 0;

  if ( (ctx->jastrow.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow.aord_num > 0);
  *aord_num = ctx->jastrow.aord_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_bord_num (const qmckl_context context, int64_t* const bord_num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  if (bord_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_bord_num",
                           "aord_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 0;

  if ( (ctx->jastrow.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow.bord_num > 0);
  *bord_num = ctx->jastrow.bord_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_cord_num (const qmckl_context context, int64_t* const cord_num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  if (cord_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_cord_num",
                           "aord_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 0;

  if ( (ctx->jastrow.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow.cord_num > 0);
  *cord_num = ctx->jastrow.cord_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_type_nucl_num (const qmckl_context context, int64_t* const type_nucl_num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  if (type_nucl_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_type_nucl_num",
                           "type_nucl_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 1;

  if ( (ctx->jastrow.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow.type_nucl_num > 0);
  *type_nucl_num = ctx->jastrow.type_nucl_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_type_nucl_vector (const qmckl_context context,
                                    int64_t* const type_nucl_vector,
                                    const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  if (type_nucl_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_type_nucl_vector",
                           "type_nucl_vector is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 2;

  if ( (ctx->jastrow.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow.type_nucl_vector != NULL);
  if (size_max < ctx->jastrow.type_nucl_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_type_nucl_vector",
                           "Array too small. Expected jastrow.type_nucl_num");
  }

  memcpy(type_nucl_vector, ctx->jastrow.type_nucl_vector, ctx->jastrow.type_nucl_num*sizeof(int64_t));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_aord_vector (const qmckl_context context,
                               double * const aord_vector,
                               const int64_t size_max) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  if (aord_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_aord_vector",
                           "aord_vector is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 3;

  if ( (ctx->jastrow.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow.aord_vector != NULL);
  int64_t sze = (ctx->jastrow.aord_num + 1)*ctx->jastrow.type_nucl_num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_aord_vector",
                           "Array too small. Expected (ctx->jastrow.aord_num + 1)*ctx->jastrow.type_nucl_num");
  }
  memcpy(aord_vector, ctx->jastrow.aord_vector, sze*sizeof(double));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_bord_vector (const qmckl_context context,
                               double * const bord_vector,
                               const int64_t size_max) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  if (bord_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_bord_vector",
                           "bord_vector is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 4;

  if ( (ctx->jastrow.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow.bord_vector != NULL);
  int64_t sze=ctx->jastrow.bord_num +1;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_bord_vector",
                           "Array too small. Expected (ctx->jastrow.bord_num + 1)");
  }
  memcpy(bord_vector, ctx->jastrow.bord_vector, sze*sizeof(double));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_cord_vector (const qmckl_context context,
                               double * const cord_vector,
                               const int64_t size_max) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  if (cord_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_cord_vector",
                           "cord_vector is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 5;

  if ( (ctx->jastrow.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow.cord_vector != NULL);

  int64_t dim_cord_vect;
  qmckl_exit_code rc = qmckl_get_jastrow_dim_cord_vect(context, &dim_cord_vect);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t sze=dim_cord_vect * ctx->jastrow.type_nucl_num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_cord_vector",
                           "Array too small. Expected dim_cord_vect * jastrow.type_nucl_num");
  }
  memcpy(cord_vector, ctx->jastrow.cord_vector, sze*sizeof(double));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_jastrow_ord_num(qmckl_context context,
                          const int64_t aord_num,
                          const int64_t bord_num,
                          const int64_t cord_num)
{

  int32_t mask = 1;

if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
  return QMCKL_NULL_CONTEXT;
 }

qmckl_context_struct* const ctx = (qmckl_context_struct*) context;


if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
    return qmckl_failwith( context,
                           QMCKL_ALREADY_SET,
                           "qmckl_set_jastrow_*",
                           NULL);

 }


  if (aord_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_ord_num",
                           "aord_num <= 0");
  }

  if (bord_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_ord_num",
                           "bord_num <= 0");
  }

  if (cord_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_ord_num",
                           "cord_num <= 0");
  }

  ctx->jastrow.aord_num = aord_num;
  ctx->jastrow.bord_num = bord_num;
  ctx->jastrow.cord_num = cord_num;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_jastrow_type_nucl_num(qmckl_context context, const int64_t type_nucl_num)
{

  int32_t mask = 1 << 1;

if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
  return QMCKL_NULL_CONTEXT;
 }

qmckl_context_struct* const ctx = (qmckl_context_struct*) context;


if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
    return qmckl_failwith( context,
                           QMCKL_ALREADY_SET,
                           "qmckl_set_jastrow_*",
                           NULL);

 }


  if (type_nucl_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_type_nucl_num",
                           "type_nucl_num < 0");
  }

  ctx->jastrow.type_nucl_num = type_nucl_num;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_jastrow_type_nucl_vector(qmckl_context context,
                                   int64_t const * type_nucl_vector,
                                   const int64_t nucl_num)
{

  int32_t mask = 1 << 2;

if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
  return QMCKL_NULL_CONTEXT;
 }

qmckl_context_struct* const ctx = (qmckl_context_struct*) context;


if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
    return qmckl_failwith( context,
                           QMCKL_ALREADY_SET,
                           "qmckl_set_jastrow_*",
                           NULL);

 }


  int64_t type_nucl_num;
  qmckl_exit_code rc = qmckl_get_jastrow_type_nucl_num(context, &type_nucl_num);
  if (rc != QMCKL_SUCCESS) return rc;

  if (type_nucl_num == 0) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_jastrow_type_nucl_vector",
                           "type_nucl_num is not set");
  }

  if (type_nucl_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_type_nucl_vector",
                           "type_nucl_vector = NULL");
  }

  if (ctx->jastrow.type_nucl_vector != NULL) {
    rc = qmckl_free(context, ctx->jastrow.type_nucl_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_type_nucl_vector",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = nucl_num * sizeof(int64_t);
  int64_t* new_array = (int64_t*) qmckl_malloc(context, mem_info);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_type_nucl_vector",
                           NULL);
  }

  memcpy(new_array, type_nucl_vector, mem_info.size);

  ctx->jastrow.type_nucl_vector = new_array;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_jastrow_aord_vector(qmckl_context context,
                              double const * aord_vector,
                              const int64_t size_max)
{
  int32_t mask = 1 << 3;

if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
  return QMCKL_NULL_CONTEXT;
 }

qmckl_context_struct* const ctx = (qmckl_context_struct*) context;


if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
    return qmckl_failwith( context,
                           QMCKL_ALREADY_SET,
                           "qmckl_set_jastrow_*",
                           NULL);

 }


  int64_t aord_num;
  qmckl_exit_code rc = qmckl_get_jastrow_aord_num(context, &aord_num);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t type_nucl_num;
  rc = qmckl_get_jastrow_type_nucl_num(context, &type_nucl_num);
  if (rc != QMCKL_SUCCESS) return rc;

  if (aord_num == 0) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_jastrow_coefficient",
                           "aord_num is not set");
  }

  if (aord_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_aord_vector",
                           "aord_vector = NULL");
  }

  if (ctx->jastrow.aord_vector != NULL) {
    rc = qmckl_free(context, ctx->jastrow.aord_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ord_vector",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = (aord_num + 1) * type_nucl_num * sizeof(double);

  if ((size_t) size_max < mem_info.size/sizeof(double)) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_aord_vector",
                           "Array too small. Expected (aord_num+1)*type_nucl_num");
  }

  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_coefficient",
                           NULL);
  }

  memcpy(new_array, aord_vector, mem_info.size);

  ctx->jastrow.aord_vector = new_array;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_jastrow_bord_vector(qmckl_context context,
                              double const * bord_vector,
                              const int64_t size_max)
{
  int32_t mask = 1 << 4;

if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
  return QMCKL_NULL_CONTEXT;
 }

qmckl_context_struct* const ctx = (qmckl_context_struct*) context;


if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
    return qmckl_failwith( context,
                           QMCKL_ALREADY_SET,
                           "qmckl_set_jastrow_*",
                           NULL);

 }


  int64_t bord_num;
  qmckl_exit_code rc = qmckl_get_jastrow_bord_num(context, &bord_num);
  if (rc != QMCKL_SUCCESS) return rc;

  if (bord_num == 0) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_jastrow_coefficient",
                           "bord_num is not set");
  }

  if (bord_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_bord_vector",
                           "bord_vector = NULL");
  }

  if (ctx->jastrow.bord_vector != NULL) {
    rc = qmckl_free(context, ctx->jastrow.bord_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ord_vector",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = (bord_num + 1) * sizeof(double);

  if ((size_t) size_max < mem_info.size/sizeof(double)) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_bord_vector",
                           "Array too small. Expected (bord_num+1)");
  }

  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_coefficient",
                           NULL);
  }

  memcpy(new_array, bord_vector, mem_info.size);

  ctx->jastrow.bord_vector = new_array;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_jastrow_cord_vector(qmckl_context context,
                              double const * cord_vector,
                              const int64_t size_max)
{
  int32_t mask = 1 << 5;

if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
  return QMCKL_NULL_CONTEXT;
 }

qmckl_context_struct* const ctx = (qmckl_context_struct*) context;


if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
    return qmckl_failwith( context,
                           QMCKL_ALREADY_SET,
                           "qmckl_set_jastrow_*",
                           NULL);

 }


  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t dim_cord_vect;
  rc = qmckl_get_jastrow_dim_cord_vect(context, &dim_cord_vect);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t type_nucl_num;
  rc = qmckl_get_jastrow_type_nucl_num(context, &type_nucl_num);
  if (rc != QMCKL_SUCCESS) return rc;

  if (dim_cord_vect == 0) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_jastrow_coefficient",
                           "dim_cord_vect is not set");
  }

  if (cord_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_cord_vector",
                           "cord_vector = NULL");
  }

  if (ctx->jastrow.cord_vector != NULL) {
    rc = qmckl_free(context, ctx->jastrow.cord_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_cord_vector",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = dim_cord_vect * type_nucl_num * sizeof(double);

  if ((size_t) size_max < mem_info.size/sizeof(double)) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_cord_vector",
                           "Array too small. Expected dim_cord_vect * type_nucl_num");
  }

  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_coefficient",
                           NULL);
  }

  memcpy(new_array, cord_vector, mem_info.size);

  ctx->jastrow.cord_vector = new_array;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_set_jastrow_type_nucl_vector_device(qmckl_context context,
                                          const int64_t * type_nucl_vector,
                                          const int64_t nucl_num,
                                          int device_id)
{

  int32_t mask = 1 << 2;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  
  if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_*",
                             NULL);
  
   }
  

  int64_t type_nucl_num;
  qmckl_exit_code rc = qmckl_get_jastrow_type_nucl_num(context, &type_nucl_num);

  if (rc != QMCKL_SUCCESS) return rc;

  if (type_nucl_num == 0) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_jastrow_type_nucl_vector",
                           "type_nucl_num is not set");
  }

  if (type_nucl_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_type_nucl_vector",
                           "type_nucl_vector = NULL");
  }

  if (ctx->jastrow.type_nucl_vector_device != NULL) {
    rc = qmckl_free(context, ctx->jastrow.type_nucl_vector_device);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_type_nucl_vector",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = nucl_num * sizeof(int64_t);
  int64_t* new_array = (int64_t*)qmckl_malloc_device(context, mem_info, device_id);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_type_nucl_vector",
                           NULL);
  }
  omp_target_memcpy(new_array, type_nucl_vector,
                    mem_info.size,
                    0, 0,
                    device_id, device_id
                    );

  ctx->jastrow.type_nucl_vector_device = new_array;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_device(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_jastrow_aord_vector_device(qmckl_context context,
                                     double const * aord_vector,
                                     const int64_t size_max,
                                     int device_id)
{
  int32_t mask = 1 << 3;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  
  if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_*",
                             NULL);
  
   }
  

  int64_t aord_num;
  qmckl_exit_code rc = qmckl_get_jastrow_aord_num(context, &aord_num);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t type_nucl_num;
  rc = qmckl_get_jastrow_type_nucl_num(context, &type_nucl_num);
  if (rc != QMCKL_SUCCESS) return rc;

  if (aord_num == 0) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_jastrow_coefficient",
                           "aord_num is not set");
  }

  if (aord_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_aord_vector",
                           "aord_vector = NULL");
  }

  if (ctx->jastrow.aord_vector != NULL) {
    rc = qmckl_free(context, ctx->jastrow.aord_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ord_vector",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = (aord_num + 1) * type_nucl_num * sizeof(double);

  if ((size_t) size_max < mem_info.size/sizeof(double)) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_aord_vector",
                           "Array too small. Expected (aord_num+1)*type_nucl_num");
  }

  double* new_array = (double*) qmckl_malloc_device(context, mem_info, device_id);


  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_coefficient",
                           NULL);
  }

  omp_target_memcpy(
    new_array, aord_vector, mem_info.size,
    0, 0,
    device_id, device_id
  );

  ctx->jastrow.aord_vector_device = new_array;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_device(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_jastrow_bord_vector_device(qmckl_context context,
                                     double const * bord_vector,
                                     const int64_t size_max,
                                     int device_id)
{
  int32_t mask = 1 << 4;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  
  if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_*",
                             NULL);
  
   }
  

  int64_t bord_num;
  qmckl_exit_code rc = qmckl_get_jastrow_bord_num(context, &bord_num);
  if (rc != QMCKL_SUCCESS) return rc;

  if (bord_num == 0) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_jastrow_coefficient",
                           "bord_num is not set");
  }

  if (bord_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_bord_vector",
                           "bord_vector = NULL");
  }

  if (ctx->jastrow.bord_vector != NULL) {
    rc = qmckl_free(context, ctx->jastrow.bord_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ord_vector",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = (bord_num + 1) * sizeof(double);

  if ((size_t) size_max < mem_info.size/sizeof(double)) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_bord_vector",
                           "Array too small. Expected (bord_num+1)");
  }

  double* new_array = (double*) qmckl_malloc_device(context, mem_info, device_id);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_coefficient",
                           NULL);
  }

  omp_target_memcpy(
    new_array, bord_vector, mem_info.size,
    0, 0,
    device_id, device_id
  );

  ctx->jastrow.bord_vector_device = new_array;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_device(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_set_jastrow_cord_vector_device(qmckl_context context,
                                     double const * cord_vector,
                                     const int64_t size_max,
                                     int device_id)
{

  int32_t mask = 1 << 5;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  
  if (mask != 0 && !(ctx->jastrow.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_*",
                             NULL);
  
   }
  


  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t dim_cord_vect;
  rc = qmckl_get_jastrow_dim_cord_vect(context, &dim_cord_vect);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t type_nucl_num;
  rc = qmckl_get_jastrow_type_nucl_num(context, &type_nucl_num);
  if (rc != QMCKL_SUCCESS) return rc;

  if (dim_cord_vect == 0) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_jastrow_coefficient",
                           "dim_cord_vect is not set");
  }

  if (cord_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_cord_vector",
                           "cord_vector = NULL");
  }

  if (ctx->jastrow.cord_vector != NULL) {
    rc = qmckl_free(context, ctx->jastrow.cord_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ord_vector",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = dim_cord_vect * type_nucl_num * sizeof(double);

  if ((size_t) size_max < mem_info.size/sizeof(double)) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_cord_vector",
                           "Array too small. Expected dim_cord_vect * type_nucl_num");
  }

  double* new_array = (double*) qmckl_malloc_device(context, mem_info, device_id);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_coefficient",
                           NULL);
  }

  omp_target_memcpy(
    new_array, cord_vector, mem_info.size,
    0, 0,
    device_id, device_id
  );

  ctx->jastrow.cord_vector_device = new_array;

  ctx->jastrow.uninitialized &= ~mask;
  ctx->jastrow.provided = (ctx->jastrow.uninitialized == 0);
  if (ctx->jastrow.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_device(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code qmckl_finalize_jastrow(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* ----------------------------------- */
  /* Check for the necessary information */
  /* ----------------------------------- */

  /* Check for the electron data
     1. elec_num
     2. ee_distances_rescaled
  */
  if (!(ctx->electron.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  /* Check for the nucleus data
     1. nucl_num
     2. en_distances_rescaled
  */
  if (!(ctx->nucleus.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_nucleus",
                           NULL);
  }

  /* Decide if the Jastrow should be offloaded on GPU or not */
#if defined(HAVE_HPC) && (defined(HAVE_CUBLAS_OFFLOAD) || defined(HAVE_OPENACC_OFFLOAD) || defined(HAVE_OPENMP_OFFLOAD))
  ctx->jastrow.gpu_offload = true; // ctx->electron.num > 100;
#endif

  qmckl_exit_code rc = QMCKL_SUCCESS;
  return rc;


}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_finalize_jastrow_device(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* ----------------------------------- */
  /* Check for the necessary information */
  /* ----------------------------------- */

  /* Check for the electron data
     1. elec_num
     2. ee_distances_rescaled
  */
  if (!(ctx->electron.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  /* Check for the nucleus data
     1. nucl_num
     2. en_distances_rescaled
  */
  if (!(ctx->nucleus.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_nucleus",
                           NULL);
  }

  qmckl_exit_code rc = QMCKL_SUCCESS;
  return rc;


}
#endif

qmckl_exit_code
qmckl_get_jastrow_asymp_jasb(qmckl_context context,
                             double* const asymp_jasb,
                             const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_asymp_jasb(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 2;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_asymp_jasb",
                           "Array too small. Expected 2");
  }
  memcpy(asymp_jasb, ctx->jastrow.asymp_jasb, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_asymp_jasb(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee kappa is provided */
  double rescale_factor_kappa_ee;
  rc = qmckl_get_electron_rescale_factor_ee(context, &rescale_factor_kappa_ee);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.asymp_jasb_date) {

    /* Allocate array */
    if (ctx->jastrow.asymp_jasb == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 2 * sizeof(double);
      double* asymp_jasb = (double*) qmckl_malloc(context, mem_info);

      if (asymp_jasb == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_asymp_jasb",
                               NULL);
      }
      ctx->jastrow.asymp_jasb = asymp_jasb;
    }

    rc = qmckl_compute_asymp_jasb(context,
                                  ctx->jastrow.bord_num,
                                  ctx->jastrow.bord_vector,
                                  rescale_factor_kappa_ee,
                                  ctx->jastrow.asymp_jasb);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.asymp_jasb_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_asymp_jasb (
	  const qmckl_context context,
	  const int64_t bord_num,
	  const double* bord_vector,
	  const double rescale_factor_kappa_ee,
	  double* const asymp_jasb ) {

    if (context == QMCKL_NULL_CONTEXT){
      return QMCKL_INVALID_CONTEXT;
    }

    if (bord_num <= 0) {
      return QMCKL_INVALID_ARG_2;
    }

    const double kappa_inv = 1.0 / rescale_factor_kappa_ee;
    const double asym_one = bord_vector[0] * kappa_inv / (1.0 + bord_vector[1] * kappa_inv);
    asymp_jasb[0] = asym_one;
    asymp_jasb[1] = 0.5 * asym_one;

    for (int i = 0 ; i <= 1; ++i) {
      double x = kappa_inv;
      for (int p = 1; p < bord_num; ++p){
        x *= kappa_inv;
        asymp_jasb[i] = asymp_jasb[i] + bord_vector[p + 1] * x;
      }
    }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_factor_ee(qmckl_context context,
                            double* const factor_ee,
                            const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_factor_ee(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze=ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_ee",
                           "Array too small. Expected walker.num");
  }
  memcpy(factor_ee, ctx->jastrow.factor_ee, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_factor_ee(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee rescaled distance is provided */
  rc = qmckl_provide_ee_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.factor_ee_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.factor_ee);
      ctx->jastrow.factor_ee = NULL;
    }
    
    /* Allocate array */
    if (ctx->jastrow.factor_ee == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* factor_ee = (double*) qmckl_malloc(context, mem_info);

      if (factor_ee == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_factor_ee",
                               NULL);
      }
      ctx->jastrow.factor_ee = factor_ee;
    }

    rc = qmckl_compute_factor_ee(context,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->electron.up_num,
                                 ctx->jastrow.bord_num,
                                 ctx->jastrow.bord_vector,
                                 ctx->electron.ee_distance_rescaled,
                                 ctx->jastrow.asymp_jasb,
                                 ctx->jastrow.factor_ee);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.factor_ee_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_factor_ee (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t up_num,
      const int64_t bord_num,
      const double* bord_vector,
      const double* ee_distance_rescaled,
      const double* asymp_jasb,
      double* const factor_ee ) {

  int ipar; // can we use a smaller integer?
  double x, x1, spin_fact, power_ser;

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  if (bord_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  for (int nw = 0; nw < walk_num; ++nw) {
    factor_ee[nw] = 0.0; // put init array here.
    for (int i = 0; i < elec_num; ++i ) {
      for (int j = 0; j < i; ++j) {
        //x = ee_distance_rescaled[j * (walk_num * elec_num) + i * (walk_num) + nw];
        x = ee_distance_rescaled[j + i * elec_num + nw*(elec_num * elec_num)];
        x1 = x;
        power_ser = 0.0;
        spin_fact = 1.0;
        ipar = 0; // index of asymp_jasb

        for (int p = 1; p < bord_num; ++p) {
          x = x * x1;
          power_ser = power_ser + bord_vector[p + 1] * x;
        }

        if(i < up_num || j >= up_num) {
          spin_fact = 0.5;
          ipar = 1;
        }

        factor_ee[nw] = factor_ee[nw] + spin_fact * bord_vector[0]  * \
                                x1 /  \
                                (1.0 + bord_vector[1] *               \
                                x1)   \
                               -asymp_jasb[ipar] + power_ser;

      }
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_factor_ee_deriv_e(qmckl_context context,
                                    double* const factor_ee_deriv_e,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_factor_ee_deriv_e(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_ee_deriv_e",
                           "Array too small. Expected 4*walk_num*elec_num");
  }

  memcpy(factor_ee_deriv_e, ctx->jastrow.factor_ee_deriv_e, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_factor_ee_deriv_e(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee rescaled distance is provided */
  rc = qmckl_provide_ee_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if ee rescaled distance deriv e is provided */
  rc = qmckl_provide_ee_distance_rescaled_deriv_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.factor_ee_deriv_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.factor_ee_deriv_e);
      ctx->jastrow.factor_ee_deriv_e = NULL;
    }
    
    /* Allocate array */
    if (ctx->jastrow.factor_ee_deriv_e == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * 4 * ctx->electron.num * sizeof(double);
      double* factor_ee_deriv_e = (double*) qmckl_malloc(context, mem_info);

      if (factor_ee_deriv_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_factor_ee_deriv_e",
                               NULL);
      }
      ctx->jastrow.factor_ee_deriv_e = factor_ee_deriv_e;
    }

    rc = qmckl_compute_factor_ee_deriv_e(context,
                                         ctx->electron.walker.num,
                                         ctx->electron.num,
                                         ctx->electron.up_num,
                                         ctx->jastrow.bord_num,
                                         ctx->jastrow.bord_vector,
                                         ctx->electron.ee_distance_rescaled,
                                         ctx->electron.ee_distance_rescaled_deriv_e,
                                         ctx->jastrow.factor_ee_deriv_e);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.factor_ee_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_factor_ee_deriv_e_hpc(
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t up_num,
      const int64_t bord_num,
      const double* bord_vector,
      const double* ee_distance_rescaled,
      const double* ee_distance_rescaled_deriv_e,
      double* const factor_ee_deriv_e ) {

  int64_t ii;
  double  pow_ser_g[3];
  double  dx[4];
  double  x, spin_fact, y;
  double  den, invden, invden2, invden3, xinv;
  double  lap1, lap2, lap3, third;

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  } 

  if (walk_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  } 

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (bord_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  
  for (int nw = 0; nw < walk_num; ++nw) {
    for (int ii = 0; ii < 4; ++ii) {
      for (int j = 0; j < elec_num; ++j) {
        factor_ee_deriv_e[j + ii * elec_num + nw * elec_num * 4]  = 0.0;
      }
    }
  }
  
  third = 1.0 / 3.0;

  for (int nw = 0; nw < walk_num; ++nw) {
    for (int i = 0; i < elec_num; ++i) {
      for (int j = 0; j < elec_num; ++j) {
        x = ee_distance_rescaled[j + i * elec_num + nw * elec_num * elec_num];
        if (fabs(x) < 1.0e-18) continue;
        for (int ii = 0; ii < 3; ++ii){
            pow_ser_g[ii] = 0.0;
        }    
        spin_fact   = 1.0;
        den         = 1.0 + bord_vector[1] * x;
        invden      = 1.0 / den;
        invden2     = invden * invden;
        invden3     = invden2 * invden;
        xinv        = 1.0 / (x + 1.0e-18);
        
        dx[0] = ee_distance_rescaled_deriv_e[0 \
                                           + j * 4 + i * 4 * elec_num \
                                           + nw * 4 * elec_num * elec_num];
        dx[1] = ee_distance_rescaled_deriv_e[1  \
                                           + j * 4 + i * 4 * elec_num \
                                           + nw * 4 * elec_num * elec_num];
        dx[2] = ee_distance_rescaled_deriv_e[2  \
                                           + j * 4 + i * 4 * elec_num \
                                           + nw * 4 * elec_num * elec_num];
        dx[3] = ee_distance_rescaled_deriv_e[3  \
                                           + j * 4 + i * 4 * elec_num \
                                           + nw * 4 * elec_num * elec_num];

        if((i <= (up_num-1) && j <= (up_num-1) ) || (i > (up_num-1) && j > (up_num-1))) {
          spin_fact = 0.5;
        }

        lap1 = 0.0;
        lap2 = 0.0;
        lap3 = 0.0;
        for (int ii = 0; ii < 3; ++ii) {
          x = ee_distance_rescaled[j + i * elec_num + nw * elec_num * elec_num];
          if (fabs(x) < 1.0e-18) continue;
          for (int p = 2; p < bord_num+1; ++p) {
            y = p * bord_vector[(p-1) + 1] * x;
            pow_ser_g[ii] = pow_ser_g[ii] + y * dx[ii];
            lap1 = lap1 + (p - 1) * y * xinv * dx[ii] * dx[ii];
            lap2 = lap2 + y;
            x = x * ee_distance_rescaled[j + i * elec_num + nw * elec_num * elec_num];
          }

          lap3 = lap3 - 2.0 * bord_vector[1] * dx[ii] * dx[ii];

          factor_ee_deriv_e[i  + ii * elec_num  + nw * elec_num * 4 ] +=              \
                                       + spin_fact * bord_vector[0] * dx[ii] * invden2 \
                                       + pow_ser_g[ii] ;
        }

        ii = 3;
        lap2 = lap2 * dx[ii] * third;
        lap3 = lap3 + den * dx[ii];
        lap3 = lap3 * (spin_fact * bord_vector[0] * invden3);
        factor_ee_deriv_e[i + ii*elec_num + nw * elec_num * 4] += lap1 + lap2 + lap3;

      }
    }
  }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_factor_ee_deriv_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t up_num,
      const int64_t bord_num,
      const double* bord_vector,
      const double* ee_distance_rescaled,
      const double* ee_distance_rescaled_deriv_e,
      double* const factor_ee_deriv_e ) {

      #ifdef HAVE_HPC
      return qmckl_compute_factor_ee_deriv_e_hpc(context, walk_num, elec_num, up_num, bord_num, bord_vector, ee_distance_rescaled, ee_distance_rescaled_deriv_e, factor_ee_deriv_e ); 
      #else
      return qmckl_compute_factor_ee_deriv_e_doc(context, walk_num, elec_num, up_num, bord_num, bord_vector, ee_distance_rescaled, ee_distance_rescaled_deriv_e, factor_ee_deriv_e ); 
      #endif
}

qmckl_exit_code
qmckl_get_jastrow_factor_en(qmckl_context context,
                            double* const factor_en,
                            const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_factor_en(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze=ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_en",
                           "Array too small. Expected walker.num");
  }
  memcpy(factor_en, ctx->jastrow.factor_en, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_factor_en(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.factor_en_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.factor_en);
      ctx->jastrow.factor_en = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.factor_en == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* factor_en = (double*) qmckl_malloc(context, mem_info);

      if (factor_en == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_factor_en",
                               NULL);
      }
      ctx->jastrow.factor_en = factor_en;
    }

    rc = qmckl_compute_factor_en(context,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->nucleus.num,
                                 ctx->jastrow.type_nucl_num,
                                 ctx->jastrow.type_nucl_vector,
                                 ctx->jastrow.aord_num,
                                 ctx->jastrow.aord_vector,
                                 ctx->electron.en_distance_rescaled,
                                 ctx->jastrow.factor_en);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }
    
    ctx->jastrow.factor_en_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_factor_en (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const int64_t aord_num,
      const double* aord_vector,
      const double* en_distance_rescaled,
      double* const factor_en ) {

  double  x, x1, power_ser;


  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  if (type_nucl_num <= 0) {
    return QMCKL_INVALID_ARG_5;
  }

  if (type_nucl_vector == NULL) {
    return QMCKL_INVALID_ARG_6;
  }

  if (aord_num <= 0) {
     return QMCKL_INVALID_ARG_7;
  }

  if (aord_vector == NULL) {
    return QMCKL_INVALID_ARG_8;
  }

  if (en_distance_rescaled == NULL) {
    return QMCKL_INVALID_ARG_9;
  }

  if (factor_en == NULL) {
    return QMCKL_INVALID_ARG_10;
  }


  for (int nw = 0; nw < walk_num; ++nw ) {
    // init array
    factor_en[nw] = 0.0;
    for (int a = 0; a < nucl_num; ++a ) {
      for (int i = 0; i < elec_num; ++i ) {
        // x = ee_distance_rescaled[j * (walk_num * elec_num) + i * (walk_num) + nw];
        x = en_distance_rescaled[i + a * elec_num + nw * (elec_num * nucl_num)];
        x1 = x;
        power_ser = 0.0;

        for (int p = 2; p < aord_num+1; ++p) {
          x = x * x1;
          power_ser = power_ser + aord_vector[(p+1)-1 + (type_nucl_vector[a]-1) * aord_num] * x;
        }

        factor_en[nw] = factor_en[nw] + aord_vector[0 + (type_nucl_vector[a]-1)*aord_num] * x1 / \
                        (1.0 + aord_vector[1 + (type_nucl_vector[a]-1) * aord_num] * x1) + \
                        power_ser;

      }
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_factor_en_deriv_e(qmckl_context context,
                                    double* const factor_en_deriv_e,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_factor_en_deriv_e(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_en_deriv_e",
                           "Array too small. Expected 4*walker.num*elec_num");
  }
  memcpy(factor_en_deriv_e, ctx->jastrow.factor_en_deriv_e, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_factor_en_deriv_e(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_en_distance_rescaled_deriv_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.factor_en_deriv_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.factor_en_deriv_e);
      ctx->jastrow.factor_en_deriv_e = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.factor_en_deriv_e == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * 4 * ctx->electron.num * sizeof(double);
      double* factor_en_deriv_e = (double*) qmckl_malloc(context, mem_info);

      if (factor_en_deriv_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_factor_en_deriv_e",
                               NULL);
      }
      ctx->jastrow.factor_en_deriv_e = factor_en_deriv_e;
    }

    rc = qmckl_compute_factor_en_deriv_e(context,
                                         ctx->electron.walker.num,
                                         ctx->electron.num,
                                         ctx->nucleus.num,
                                         ctx->jastrow.type_nucl_num,
                                         ctx->jastrow.type_nucl_vector,
                                         ctx->jastrow.aord_num,
                                         ctx->jastrow.aord_vector,
                                         ctx->electron.en_distance_rescaled,
                                         ctx->electron.en_distance_rescaled_deriv_e,
                                         ctx->jastrow.factor_en_deriv_e);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.factor_en_deriv_e_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_een_rescaled_e(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_e(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.num * ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_een_rescaled_e",
                           "Array too small. Expected ctx->electron.num * ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow.cord_num + 1)");
  }
  memcpy(distance_rescaled, ctx->jastrow.een_rescaled_e, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_e(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_ee_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.een_rescaled_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.een_rescaled_e);
      ctx->jastrow.een_rescaled_e = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.een_rescaled_e == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * ctx->electron.num *
        ctx->electron.walker.num * (ctx->jastrow.cord_num + 1) * sizeof(double);
      double* een_rescaled_e = (double*) qmckl_malloc(context, mem_info);

      if (een_rescaled_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_een_rescaled_e",
                               NULL);
      }
      ctx->jastrow.een_rescaled_e = een_rescaled_e;
    }

    rc = qmckl_compute_een_rescaled_e(context,
                                      ctx->electron.walker.num,
                                      ctx->electron.num,
                                      ctx->jastrow.cord_num,
                                      ctx->electron.rescale_factor_kappa_ee,
                                      ctx->electron.ee_distance,
                                      ctx->jastrow.een_rescaled_e);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.een_rescaled_e_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_een_rescaled_e_device(qmckl_context context, int device_id)
{


  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_ee_distance_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.een_rescaled_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      omp_target_free(ctx->jastrow.een_rescaled_e_device, device_id);
      ctx->jastrow.een_rescaled_e_device = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.een_rescaled_e_device == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * ctx->electron.num *
        ctx->electron.walker.num * (ctx->jastrow.cord_num + 1) * sizeof(double);
      double* een_rescaled_e = (double*) qmckl_malloc_device(context, mem_info, device_id);

      if (een_rescaled_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_een_rescaled_e_device",
                               NULL);
      }
      ctx->jastrow.een_rescaled_e_device = een_rescaled_e;
    }

    rc = qmckl_compute_een_rescaled_e_device(context,
                                      ctx->electron.walker.num,
                                      ctx->electron.num,
                                      ctx->jastrow.cord_num,
                                      ctx->electron.rescale_factor_kappa_ee,
                                      ctx->electron.ee_distance_device,
                                      ctx->jastrow.een_rescaled_e_device,
                                      device_id);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.een_rescaled_e_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code qmckl_compute_een_rescaled_e_hpc (
	  const qmckl_context context,
	  const int64_t walk_num,
	  const int64_t elec_num,
	  const int64_t cord_num,
	  const double rescale_factor_kappa_ee,
	  const double* ee_distance,
	  double* const een_rescaled_e ) {

  double   *een_rescaled_e_ij;
  double   x;
  const int64_t   elec_pairs = (elec_num * (elec_num - 1)) / 2;
  const int64_t   len_een_ij = elec_pairs * (cord_num + 1);
  int64_t   k;

  // number of element for the een_rescaled_e_ij[N_e*(N_e-1)/2][cord+1]
  // probably in C is better [cord+1, Ne*(Ne-1)/2]
  //elec_pairs = (elec_num * (elec_num - 1)) / 2;
  //len_een_ij = elec_pairs * (cord_num + 1);
  een_rescaled_e_ij = (double *) malloc (len_een_ij * sizeof(double));

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  if (cord_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  // Prepare table of exponentiated distances raised to appropriate power
  // init

  for (int kk = 0; kk < walk_num*(cord_num+1)*elec_num*elec_num; ++kk) {
    een_rescaled_e[kk]= 0.0;
  }

  /*
  for (int nw = 0; nw < walk_num; ++nw) {
    for (int l = 0; l < (cord_num + 1); ++l) {
      for (int i = 0; i < elec_num; ++i) {
        for (int j = 0; j < elec_num; ++j) {
          een_rescaled_e[j + i*elec_num + l*elec_num*elec_num + nw*(cord_num+1)*elec_num*elec_num]= 0.0;
        }
      }
    }
  }
  */

  for (int nw = 0; nw < walk_num; ++nw) {

    for (int kk = 0; kk < len_een_ij; ++kk) {
      // this array initialized at 0 except een_rescaled_e_ij(:, 1) = 1.0d0
      // and the arrangement of indices is [cord_num+1, ne*(ne-1)/2]
      een_rescaled_e_ij[kk]= ( kk < (elec_pairs) ? 1.0 : 0.0 );
    }

    k = 0;
    for (int i = 0; i < elec_num; ++i) {
      for (int j = 0; j < i; ++j) {
        // een_rescaled_e_ij(k, 2) = dexp(-rescale_factor_kappa_ee * ee_distance(i, j, nw));
        een_rescaled_e_ij[k + elec_pairs] = exp(-rescale_factor_kappa_ee * \
                                    ee_distance[j + i*elec_num + nw*(elec_num*elec_num)]);
        k = k + 1;
      }
    }


    for (int l = 2; l < (cord_num+1); ++l) {
      for (int k = 0; k < elec_pairs; ++k) {
      // een_rescaled_e_ij(k, l + 1) = een_rescaled_e_ij(k, l + 1 - 1) * een_rescaled_e_ij(k, 2)
        een_rescaled_e_ij[k+l*elec_pairs] = een_rescaled_e_ij[k + (l - 1)*elec_pairs] * \
                                                  een_rescaled_e_ij[k + elec_pairs];
      }
    }


  // prepare the actual een table
  for (int i = 0; i < elec_num; ++i){
    for (int j = 0; j < elec_num; ++j) {
      een_rescaled_e[j + i*elec_num + 0 + nw*(cord_num+1)*elec_num*elec_num] = 1.0;
    }
  }

    // Up to here it should work.
  for ( int l = 1; l < (cord_num+1); ++l) {
    k = 0;
    for (int i = 0; i < elec_num; ++i) {
      for (int j = 0; j < i; ++j) {
        x = een_rescaled_e_ij[k + l*elec_pairs];
        een_rescaled_e[j + i*elec_num + l*elec_num*elec_num + nw*elec_num*elec_num*(cord_num+1)] = x;
        een_rescaled_e[i + j*elec_num + l*elec_num*elec_num + nw*elec_num*elec_num*(cord_num+1)] = x;
        k = k + 1;
      }
    }
  }

  for (int l = 0; l < (cord_num + 1); ++l) {
    for (int j = 0; j < elec_num; ++j) {
      een_rescaled_e[j + j*elec_num + l*elec_num*elec_num + nw*elec_num*elec_num*(cord_num+1)] = 0.0;
    }
  }

  }

  free(een_rescaled_e_ij);

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_een_rescaled_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_ee,
      const double* ee_distance,
      double* const een_rescaled_e ) {

#ifdef HAVE_HPC
return qmckl_compute_een_rescaled_e_hpc(context, walk_num, elec_num, cord_num, rescale_factor_kappa_ee, ee_distance, een_rescaled_e);
#else
return qmckl_compute_een_rescaled_e_doc(context, walk_num, elec_num, cord_num, rescale_factor_kappa_ee, ee_distance, een_rescaled_e);
#endif
}

/* Device pointers */


qmckl_exit_code qmckl_compute_een_rescaled_e_device (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t cord_num,
  const double rescale_factor_kappa_ee,
  const double* ee_distance,
  double* const een_rescaled_e,
  int device_id ) {


  double   *een_rescaled_e_ij;
  double   x;
  const int64_t   elec_pairs = (elec_num * (elec_num - 1)) / 2;
  const int64_t   len_een_ij = elec_pairs * (cord_num + 1);
  int64_t   k;

  // number of element for the een_rescaled_e_ij[N_e*(N_e-1)/2][cord+1]
  // probably in C is better [cord+1, Ne*(Ne-1)/2]
  //elec_pairs = (elec_num * (elec_num - 1)) / 2;
  //len_een_ij = elec_pairs * (cord_num + 1);
  een_rescaled_e_ij = (double *) omp_target_alloc (len_een_ij * sizeof(double), device_id);

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  if (cord_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  // Prepare table of exponentiated distances raised to appropriate power
  // init

  // TODO Parallelism can probably be heavily optimized

  #pragma omp target is_device_ptr(een_rescaled_e)
  {
  #pragma omp teams distribute parallel for simd
  for (int kk = 0; kk < walk_num*(cord_num+1)*elec_num*elec_num; ++kk) {
    een_rescaled_e[kk]= 0.0;
  }
  }

  #pragma omp target is_device_ptr(ee_distance, een_rescaled_e, een_rescaled_e_ij)
  {
  for (int nw = 0; nw < walk_num; ++nw) {
    #pragma omp parallel for
    for (int kk = 0; kk < len_een_ij; ++kk) {
      // this array initialized at 0 except een_rescaled_e_ij(:, 1) = 1.0d0
      // and the arrangement of indices is [cord_num+1, ne*(ne-1)/2]
      een_rescaled_e_ij[kk]= ( kk < (elec_pairs) ? 1.0 : 0.0 );
    }

    k = 0;
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < elec_num; ++i) {
      for (int j = 0; j < i; ++j) {
        // een_rescaled_e_ij(k, 2) = dexp(-rescale_factor_kappa_ee * ee_distance(i, j, nw));
        een_rescaled_e_ij[k + elec_pairs] = exp(-rescale_factor_kappa_ee * \
                                    ee_distance[j + i*elec_num + nw*(elec_num*elec_num)]);
        k = k + 1;
      }
    }

    #pragma omp parallel for collapse(2)
    for (int l = 2; l < (cord_num+1); ++l) {
      for (int k = 0; k < elec_pairs; ++k) {
      // een_rescaled_e_ij(k, l + 1) = een_rescaled_e_ij(k, l + 1 - 1) * een_rescaled_e_ij(k, 2)
        een_rescaled_e_ij[k+l*elec_pairs] = een_rescaled_e_ij[k + (l - 1)*elec_pairs] * \
                                                  een_rescaled_e_ij[k + elec_pairs];
      }
    }

  #pragma omp parallel for collapse(2)
  // prepare the actual een table
  for (int i = 0; i < elec_num; ++i){
    for (int j = 0; j < elec_num; ++j) {
      een_rescaled_e[j + i*elec_num + 0 + nw*(cord_num+1)*elec_num*elec_num] = 1.0;
    }
  }

    // Up to here it should work.
  #pragma omp parallel for collapse(3)
  for ( int l = 1; l < (cord_num+1); ++l) {
    //k = 0;
    for (int i = 0; i < elec_num; ++i) {
      for (int j = 0; j < i; ++j) {
        x = een_rescaled_e_ij[i*elec_num+j + l*elec_pairs];
        een_rescaled_e[j + i*elec_num + l*elec_num*elec_num + nw*elec_num*elec_num*(cord_num+1)] = x;
        een_rescaled_e[i + j*elec_num + l*elec_num*elec_num + nw*elec_num*elec_num*(cord_num+1)] = x;
        //k = k + 1;
      }
    }
  }

  #pragma omp parallel for collapse(2)
  for (int l = 0; l < (cord_num + 1); ++l) {
    for (int j = 0; j < elec_num; ++j) {
      een_rescaled_e[j + j*elec_num + l*elec_num*elec_num + nw*elec_num*elec_num*(cord_num+1)] = 0.0;
    }
  }

  }
  }

  omp_target_free(een_rescaled_e_ij, device_id);

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_een_rescaled_e_deriv_e(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_e_deriv_e(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.num * 4 * ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_een_deriv_e",
                           "Array too small. Expected ctx->electron.num * 4 * ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow.cord_num + 1)");
  }
  memcpy(distance_rescaled, ctx->jastrow.een_rescaled_e_deriv_e, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_e_deriv_e(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_een_rescaled_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.een_rescaled_e_deriv_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.een_rescaled_e_deriv_e);
      ctx->jastrow.een_rescaled_e_deriv_e = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.een_rescaled_e_deriv_e == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * 4 * ctx->electron.num *
        ctx->electron.walker.num * (ctx->jastrow.cord_num + 1) * sizeof(double);
      double* een_rescaled_e_deriv_e = (double*) qmckl_malloc(context, mem_info);

      if (een_rescaled_e_deriv_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_een_rescaled_e_deriv_e",
                               NULL);
      }
      ctx->jastrow.een_rescaled_e_deriv_e = een_rescaled_e_deriv_e;
    }

    rc = qmckl_compute_factor_een_rescaled_e_deriv_e(context,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->jastrow.cord_num,
                                                     ctx->electron.rescale_factor_kappa_ee,
                                                     ctx->electron.walker.point.coord.data,
                                                     ctx->electron.ee_distance,
                                                     ctx->jastrow.een_rescaled_e,
                                                     ctx->jastrow.een_rescaled_e_deriv_e);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.een_rescaled_e_deriv_e_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_een_rescaled_e_deriv_e_device(qmckl_context context, int device_id)
{


  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_een_rescaled_e_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.een_rescaled_e_deriv_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      omp_target_free(ctx->jastrow.een_rescaled_e_deriv_e_device, device_id);
      ctx->jastrow.een_rescaled_e_deriv_e_device = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.een_rescaled_e_deriv_e_device == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * 4 * ctx->electron.num *
        ctx->electron.walker.num * (ctx->jastrow.cord_num + 1) * sizeof(double);
      double* een_rescaled_e_deriv_e = (double*) qmckl_malloc_device(context, mem_info, device_id);

      if (een_rescaled_e_deriv_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_een_rescaled_e_deriv_e_device",
                               NULL);
      }
      ctx->jastrow.een_rescaled_e_deriv_e_device = een_rescaled_e_deriv_e;
    }

    rc = qmckl_compute_factor_een_rescaled_e_deriv_e_device(context,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->jastrow.cord_num,
                                                     ctx->electron.rescale_factor_kappa_ee,
                                                     ctx->electron.walker.point.coord.data_device,
                                                     ctx->electron.ee_distance_device,
                                                     ctx->jastrow.een_rescaled_e_device,
                                                     ctx->jastrow.een_rescaled_e_deriv_e_device,
                                                     device_id);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.een_rescaled_e_deriv_e_date = ctx->date;
  }


  return QMCKL_SUCCESS;
}
#endif

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_compute_factor_een_rescaled_e_deriv_e_device (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t cord_num,
  const double rescale_factor_kappa_ee,
  const double* coord_new,
  const double* ee_distance,
  const double* een_rescaled_e,
  double* const een_rescaled_e_deriv_e,
  int device_id
){


  double rij_inv, kappa_l;

  double* elec_dist_deriv_e = omp_target_alloc(sizeof(double*) * elec_num * elec_num * 4, device_id);


  if (context == QMCKL_NULL_CONTEXT){
    return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0){
    return QMCKL_INVALID_ARG_2;
  }


  if(elec_num <= 0){
    return QMCKL_INVALID_ARG_3;
  }

  if(cord_num <= 0){
    return QMCKL_INVALID_ARG_4;
  }


  for(int nw = 0; nw < walk_num; nw++){

  #pragma omp target is_device_ptr(elec_dist_deriv_e, coord_new, ee_distance, een_rescaled_e, een_rescaled_e_deriv_e)
  {

    #pragma omp teams distribute parallel for simd
    for(int j = 0; j < elec_num; j++){
      for(int i = 0; i < elec_num; i++){
        rij_inv = 1/ee_distance[nw * elec_num * elec_num + j *elec_num + i];
        for(int ii=0; ii < 3; ii++){
          elec_dist_deriv_e[j*elec_num * 4 + i * 4 + ii] = (coord_new[nw * 3 *elec_num + ii * elec_num + i] - coord_new[nw * 3 *elec_num + ii * elec_num + j]) * rij_inv;

        }
        elec_dist_deriv_e[j*elec_num*4 + i*4 + 3] = 2 * rij_inv;
      }
      elec_dist_deriv_e[j*elec_num*4 + j*4 + 0] = 0;
      elec_dist_deriv_e[j*elec_num*4 + j*4 + 1] = 0;
      elec_dist_deriv_e[j*elec_num*4 + j*4 + 2] = 0;
      elec_dist_deriv_e[j*elec_num*4 + j*4 + 3] = 0;
    }

  }

  #pragma omp target is_device_ptr(elec_dist_deriv_e, coord_new, ee_distance, een_rescaled_e, een_rescaled_e_deriv_e)
  {

    #pragma omp teams distribute parallel for simd collapse (3)
    for(int l=1; l < cord_num+1; l++){
      for(int j=0; j < elec_num; j++){
        for(int i=0; i < elec_num ; i++){

          kappa_l = - l * rescale_factor_kappa_ee;

          een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 0 * elec_num + i] = kappa_l * elec_dist_deriv_e[j * elec_num * 4 + i * 4 + 0];


          een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 1 * elec_num + i] = \
            kappa_l * elec_dist_deriv_e[j * elec_num * 4 + i * 4 + 1];

          een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 2 * elec_num + i] = \
            kappa_l * elec_dist_deriv_e[j * elec_num * 4 + i * 4 + 2];

          een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 3 * elec_num + i] = \
            kappa_l * elec_dist_deriv_e[j * elec_num * 4 + i * 4 + 3];



          een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 3 * elec_num + i] = \
            een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 3 * elec_num + i] \
             + een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 0 * elec_num + i] * een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 0 * elec_num + i] \
             + een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 1 * elec_num + i] * een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 1 * elec_num + i] \
             + een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 2 * elec_num + i] * een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 2 * elec_num + i];


          een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 0 * elec_num + i] = een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 0 * elec_num + i] * een_rescaled_e[nw * ((cord_num+1)) * elec_num * elec_num + l * elec_num * elec_num + j * elec_num + i];

          een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 1 * elec_num + i] = een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 1 * elec_num + i] * een_rescaled_e[nw * ((cord_num+1)) * elec_num * elec_num + l * elec_num * elec_num + j * elec_num + i];

          een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 2 * elec_num + i] = een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 2 * elec_num + i] * een_rescaled_e[nw * ((cord_num+1)) * elec_num * elec_num + l * elec_num * elec_num + j * elec_num + i];

          een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 3 * elec_num + i] = een_rescaled_e_deriv_e[nw * (cord_num+1) * elec_num * 4 *elec_num + l * elec_num * 4 * elec_num + j * 4 * elec_num + 3 * elec_num + i] * een_rescaled_e[nw * ((cord_num+1)) * elec_num * elec_num + l * elec_num * elec_num + j * elec_num + i];
        }
      }
    }
  }

  }

  omp_target_free(elec_dist_deriv_e, device_id);

  return QMCKL_SUCCESS;

}
#endif

qmckl_exit_code
qmckl_get_jastrow_een_rescaled_n(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_n(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_een_deriv_e",
                           "Array too small. Expected ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow.cord_num + 1)");
  }
  memcpy(distance_rescaled, ctx->jastrow.een_rescaled_n, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_n(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.een_rescaled_n_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.een_rescaled_n);
      ctx->jastrow.een_rescaled_n = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.een_rescaled_n == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * ctx->nucleus.num *
        ctx->electron.walker.num * (ctx->jastrow.cord_num + 1) * sizeof(double);
      double* een_rescaled_n = (double*) qmckl_malloc(context, mem_info);

      if (een_rescaled_n == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_een_rescaled_n",
                               NULL);
      }
      ctx->jastrow.een_rescaled_n = een_rescaled_n;
    }

    rc = qmckl_compute_een_rescaled_n(context,
                                      ctx->electron.walker.num,
                                      ctx->electron.num,
                                      ctx->nucleus.num,
                                      ctx->jastrow.cord_num,
                                      ctx->electron.rescale_factor_kappa_en,
                                      ctx->electron.en_distance,
                                      ctx->jastrow.een_rescaled_n);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.een_rescaled_n_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_een_rescaled_n_device(qmckl_context context, int device_id)
{


  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_en_distance_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.een_rescaled_n_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      omp_target_free(ctx->jastrow.een_rescaled_n_device, device_id);
      ctx->jastrow.een_rescaled_n_device = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.een_rescaled_n_device == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * ctx->nucleus.num *
        ctx->electron.walker.num * (ctx->jastrow.cord_num + 1) * sizeof(double);
      double* een_rescaled_n = (double*) qmckl_malloc_device(context, mem_info, device_id);

      if (een_rescaled_n == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_een_rescaled_n_device",
                               NULL);
      }
      ctx->jastrow.een_rescaled_n_device = een_rescaled_n;
    }

    rc = qmckl_compute_een_rescaled_n_device(context,
                                      ctx->electron.walker.num,
                                      ctx->electron.num,
                                      ctx->nucleus.num,
                                      ctx->jastrow.cord_num,
                                      ctx->electron.rescale_factor_kappa_en,
                                      ctx->electron.en_distance_device,
                                      ctx->jastrow.een_rescaled_n_device,
                                      device_id);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.een_rescaled_n_date = ctx->date;
  }


  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code qmckl_compute_een_rescaled_n (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_en,
      const double* en_distance,
      double* const een_rescaled_n ) {


  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  if (cord_num <= 0) {
    return QMCKL_INVALID_ARG_5;
  }

  // Prepare table of exponentiated distances raised to appropriate power
  for (int i = 0; i < (walk_num*(cord_num+1)*nucl_num*elec_num); ++i) {
    een_rescaled_n[i] = 17.0;
  }

  for (int nw = 0; nw < walk_num; ++nw) {
    for (int a = 0; a < nucl_num; ++a) {
      for (int i = 0; i < elec_num; ++i) {
        // prepare the actual een table
        //een_rescaled_n(:, :, 0, nw) = 1.0d0
        een_rescaled_n[i + a * elec_num + 0 + nw * elec_num*nucl_num*(cord_num+1)] = 1.0;
        //een_rescaled_n(i, a, 1, nw) = dexp(-rescale_factor_kappa_en * en_distance(i, a, nw))
        een_rescaled_n[i + a*elec_num + elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)] = exp(-rescale_factor_kappa_en * \
                                                                              en_distance[i + a*elec_num + nw*elec_num*nucl_num]);
      }
    }

    for (int l = 2; l < (cord_num+1); ++l){
      for (int a = 0; a < nucl_num; ++a) {
        for (int i = 0; i < elec_num; ++i) {
          een_rescaled_n[i + a*elec_num + l*elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)] = een_rescaled_n[i + a*elec_num + (l-1)*elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)] * een_rescaled_n[i+a*elec_num+elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)];
        }
      }
    }

  }

  return QMCKL_SUCCESS;
}

/* Device pointers */


qmckl_exit_code qmckl_compute_een_rescaled_n_device (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const double rescale_factor_kappa_en,
      const double* en_distance,
      double* const een_rescaled_n,
      int device_id ) {


  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  if (cord_num <= 0) {
    return QMCKL_INVALID_ARG_5;
  }


  // TODO Parallelize this stuff

  #pragma omp target is_device_ptr(en_distance, een_rescaled_n)
  {
  // Prepare table of exponentiated distances raised to appropriate power
  #pragma omp teams distribute parallel for simd
  for (int i = 0; i < (walk_num*(cord_num+1)*nucl_num*elec_num); ++i) {
    een_rescaled_n[i] = 17.0;
  }
  }


  #pragma omp target is_device_ptr(en_distance, een_rescaled_n)
  {

  #pragma omp teams distribute parallel for simd
  for (int nw = 0; nw < walk_num; ++nw) {

    for (int a = 0; a < nucl_num; ++a) {
      for (int i = 0; i < elec_num; ++i) {
        // prepare the actual een table
        //een_rescaled_n(:, :, 0, nw) = 1.0d0
        een_rescaled_n[i + a * elec_num + nw * elec_num*nucl_num*(cord_num+1)] = 1.0;
        //een_rescaled_n(i, a, 1, nw) = dexp(-rescale_factor_kappa_en * en_distance(i, a, nw))
        een_rescaled_n[i + a*elec_num + elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)] = exp(-rescale_factor_kappa_en * \
                                                                              en_distance[i + a*elec_num + nw*elec_num*nucl_num]);
      }
    }

    for (int l = 2; l < (cord_num+1); ++l){
      for (int a = 0; a < nucl_num; ++a) {
        for (int i = 0; i < elec_num; ++i) {
          een_rescaled_n[i + a*elec_num + l*elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)] = een_rescaled_n[i + a*elec_num + (l-1)*elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)] * een_rescaled_n[i+a*elec_num+elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)];
        }
      }
    }

  }

  }


  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_een_rescaled_n_deriv_e(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_n_deriv_e(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.num * 4 * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_een_deriv_e",
                           "Array too small. Expected ctx->electron.num * 4 * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow.cord_num + 1)");
  }
  memcpy(distance_rescaled, ctx->jastrow.een_rescaled_n_deriv_e, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_n_deriv_e(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if ee distance is provided */
  rc = qmckl_provide_een_rescaled_n(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.een_rescaled_n_deriv_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.een_rescaled_n_deriv_e);
      ctx->jastrow.een_rescaled_n_deriv_e = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.een_rescaled_n_deriv_e == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * 4 * ctx->nucleus.num *
        ctx->electron.walker.num * (ctx->jastrow.cord_num + 1) * sizeof(double);
      double* een_rescaled_n_deriv_e = (double*) qmckl_malloc(context, mem_info);

      if (een_rescaled_n_deriv_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_een_rescaled_n_deriv_e",
                               NULL);
      }
      ctx->jastrow.een_rescaled_n_deriv_e = een_rescaled_n_deriv_e;
    }

    rc = qmckl_compute_factor_een_rescaled_n_deriv_e(context,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->nucleus.num,
                                                     ctx->jastrow.cord_num,
                                                     ctx->electron.rescale_factor_kappa_en,
                                                     ctx->electron.walker.point.coord.data,
                                                     ctx->nucleus.coord.data,
                                                     ctx->electron.en_distance,
                                                     ctx->jastrow.een_rescaled_n,
                                                     ctx->jastrow.een_rescaled_n_deriv_e);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.een_rescaled_n_deriv_e_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_een_rescaled_n_deriv_e_device(qmckl_context context, int device_id)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_en_distance_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if ee distance is provided */
  rc = qmckl_provide_een_rescaled_n_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.een_rescaled_n_deriv_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      omp_target_free(ctx->jastrow.een_rescaled_n_deriv_e_device, device_id);
      ctx->jastrow.een_rescaled_n_deriv_e_device = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.een_rescaled_n_deriv_e_device == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * 4 * ctx->nucleus.num *
        ctx->electron.walker.num * (ctx->jastrow.cord_num + 1) * sizeof(double);
      double* een_rescaled_n_deriv_e = (double*) qmckl_malloc_device(context, mem_info, device_id);

      if (een_rescaled_n_deriv_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_een_rescaled_n_deriv_e_device",
                               NULL);
      }
      ctx->jastrow.een_rescaled_n_deriv_e_device = een_rescaled_n_deriv_e;
    }

    rc = qmckl_compute_factor_een_rescaled_n_deriv_e_device(context,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->nucleus.num,
                                                     ctx->jastrow.cord_num,
                                                     ctx->electron.rescale_factor_kappa_en,
                                                     ctx->electron.walker.point.coord.data_device,
                                                     ctx->nucleus.coord.data_device,
                                                     ctx->electron.en_distance_device,
                                                     ctx->jastrow.een_rescaled_n_device,
                                                     ctx->jastrow.een_rescaled_n_deriv_e_device,
                                                     device_id);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.een_rescaled_n_deriv_e_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}
#endif

#ifdef HAVE_DEVICE_POINTERS
    qmckl_exit_code qmckl_compute_factor_een_rescaled_n_deriv_e_device (
          qmckl_context context,
          int64_t walk_num,
          int64_t elec_num,
          int64_t nucl_num,
          int64_t cord_num,
          double rescale_factor_kappa_en,
          double* coord_ee,
          double* coord_en,
          double* en_distance,
          double* een_rescaled_n,
          double* een_rescaled_n_deriv_e,
          int device_id )
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  double * elnuc_dist_deriv_e = (double*) omp_target_alloc(4*elec_num*nucl_num, device_id);

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  if (cord_num <= 0) {
     return QMCKL_INVALID_ARG_5;
  }

  /*
  double precision      , intent(in)  :: coord_ee(elec_num,3,walk_num)
  double precision      , intent(in)  :: coord_en(nucl_num,3)
  double precision      , intent(in)  :: en_distance(elec_num,nucl_num,walk_num)
  double precision      , intent(in)  :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  double precision      , intent(out) :: een_rescaled_n_deriv_e(elec_num,4,nucl_num,0:cord_num,walk_num)
  double precision,dimension(4,elec_num,nucl_num),allocatable :: elnuc_dist_deriv_e
  */

  // Pre-compute offsets for large dimensions arrays
  // For een_rescaled_n :
  int r_offset_1 = elec_num;
  int r_offset_2 = r_offset_1 * nucl_num;
  int r_offset_3 = r_offset_2 * (cord_num+1);
  // For een_rescaled_n_deriv_e:
  int rd_offset_1 = elec_num;
  int rd_offset_2 = rd_offset_1 * 4;
  int rd_offset_3 = rd_offset_2 * nucl_num;
  int rd_offset_4 = rd_offset_3 * (cord_num+1);

  // TODO Improve parallelism ?

  #pragma omp target is_device_ptr(coord_ee, coord_en, en_distance, een_rescaled_n, een_rescaled_n_deriv_e, elnuc_dist_deriv_e)
  {
  // Prepare table of exponentiated distances raised to appropriate power
  #pragma omp teams distribute parallel for simd
  for(int i=0; i<elec_num*4*nucl_num*(cord_num+1)*walk_num; i++) {
    een_rescaled_n_deriv_e[i] = 0.0;
  }
  }

  #pragma omp target is_device_ptr(coord_ee, coord_en, en_distance, een_rescaled_n, een_rescaled_n_deriv_e, elnuc_dist_deriv_e)
  {

  for(int nw=0; nw<walk_num; nw++) {

  // prepare the actual een table
  #pragma omp parallel for collapse (2)
  for(int a=0; a<nucl_num; a++) {
    for(int i=0; i<elec_num; i++) {
      double ria_inv = 1.0 / en_distance[i + a+ nw];
      for(int ii=0; ii<3; ii++) {
        elnuc_dist_deriv_e[ii + i*(4) + a*(4*elec_num)] = (coord_ee[i + ii*(elec_num) + nw*(elec_num*3)] - coord_en[a + ii*(nucl_num)]) * ria_inv;
      }
      elnuc_dist_deriv_e[3 + i*(4) + a*(4*elec_num)] = 2.0 * ria_inv;
    }
  }

  #pragma omp parallel for collapse(3)
  for(int l=0; l<(cord_num+1); l++) {
    for(int a=0; a<nucl_num; a++) {
      for(int i=0; i<elec_num; i++) {

        double kappa_l = - ((double) l) * rescale_factor_kappa_en;

        een_rescaled_n_deriv_e[i + 0*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] = kappa_l * elnuc_dist_deriv_e[0 + i*4 + a*4*elec_num];
        een_rescaled_n_deriv_e[i + 1*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] = kappa_l * elnuc_dist_deriv_e[1 + i*4 + a*4*elec_num];
        een_rescaled_n_deriv_e[i + 2*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] = kappa_l * elnuc_dist_deriv_e[2 + i*4 + a*4*elec_num];
        een_rescaled_n_deriv_e[i + 3*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] = kappa_l * elnuc_dist_deriv_e[3 + i*4 + a*4*elec_num];

        een_rescaled_n_deriv_e[i + 3*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] = een_rescaled_n_deriv_e[i + 3*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] + een_rescaled_n_deriv_e[i + 0*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] * een_rescaled_n_deriv_e[i + 0*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] + een_rescaled_n_deriv_e[i + 1*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] * een_rescaled_n_deriv_e[i + 1*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] + een_rescaled_n_deriv_e[i + 2*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] * een_rescaled_n_deriv_e[i + 2*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3+ nw*rd_offset_4];

        een_rescaled_n_deriv_e[i + 0*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] = een_rescaled_n_deriv_e[i + 0*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] * een_rescaled_n[i + a*r_offset_1 + l*r_offset_2 + nw*r_offset_3];
        een_rescaled_n_deriv_e[i + 1*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] = een_rescaled_n_deriv_e[i + 1*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_3] * een_rescaled_n[i + a*r_offset_1 + l*r_offset_2 + nw*r_offset_3];
        een_rescaled_n_deriv_e[i + 2*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] = een_rescaled_n_deriv_e[i + 2*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] * een_rescaled_n[i + a*r_offset_1 + l*r_offset_1 + nw*r_offset_3];
        een_rescaled_n_deriv_e[i + 3*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] = een_rescaled_n_deriv_e[i + 3*rd_offset_1 + a*rd_offset_2 + l*rd_offset_3 + nw*rd_offset_4] * een_rescaled_n[i + a*r_offset_1 + l*r_offset_2 + nw*r_offset_3];
      }
    }
  }
  }

  }

  omp_target_free(elnuc_dist_deriv_e, device_id);
  return rc;
}
#endif

qmckl_exit_code qmckl_get_jastrow_dim_cord_vect(qmckl_context context, int64_t* const dim_cord_vect)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  *dim_cord_vect = ctx->jastrow.dim_cord_vect;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_cord_vect_full(qmckl_context context, double* const cord_vect_full)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_cord_vect_full(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->jastrow.dim_cord_vect * ctx->nucleus.num;
  memcpy(cord_vect_full, ctx->jastrow.cord_vect_full, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_lkpm_combined_index(qmckl_context context, int64_t* const lkpm_combined_index)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_cord_vect_full(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->jastrow.dim_cord_vect * 4;
  memcpy(lkpm_combined_index, ctx->jastrow.lkpm_combined_index, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_tmp_c(qmckl_context context, double* const tmp_c)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_cord_vect_full(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_tmp_c(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1)
               * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;
  memcpy(tmp_c, ctx->jastrow.tmp_c, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_dtmp_c(qmckl_context context, double* const dtmp_c)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_cord_vect_full(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_dtmp_c(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1)
               *4* ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;
  memcpy(dtmp_c, ctx->jastrow.dtmp_c, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_cord_vect_full_device(qmckl_context context, double* const cord_vect_full, int device_id)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_cord_vect_full_device(context, device_id);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->jastrow.dim_cord_vect * ctx->nucleus.num;
  // memcpy(cord_vect_full, ctx->jastrow.cord_vect_full, sze * sizeof(double));
  omp_target_memcpy(
    cord_vect_full, ctx->jastrow.cord_vect_full_device,
    sze * sizeof(double),
    0, 0,
    device_id, device_id
  );

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_lkpm_combined_index_device(qmckl_context context, int64_t* const lkpm_combined_index, int device_id)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_cord_vect_full_device(context, device_id);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->jastrow.dim_cord_vect * 4;
  // memcpy(lkpm_combined_index, ctx->jastrow.lkpm_combined_index, sze * sizeof(double));
  omp_target_memcpy(
    lkpm_combined_index, ctx->jastrow.lkpm_combined_index,
    sze * sizeof(double),
    0, 0,
    device_id, device_id
  );

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_tmp_c_device(qmckl_context context, double* const tmp_c, int device_id)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_cord_vect_full_device(context, device_id);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_tmp_c_device(context, device_id);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1)
               * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;
  // memcpy(tmp_c, ctx->jastrow.tmp_c, sze * sizeof(double));
  omp_target_memcpy(
    tmp_c, ctx->jastrow.tmp_c_device,
    sze * sizeof(double),
    0, 0,
    device_id, device_id
  );


  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_dtmp_c_device(qmckl_context context, double* const dtmp_c, int device_id)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_dim_cord_vect(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_cord_vect_full_device(context, device_id);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_dtmp_c_device(context, device_id);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1)
               *4* ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;
  // memcpy(dtmp_c, ctx->jastrow.dtmp_c, sze * sizeof(double));
  omp_target_memcpy(
    dtmp_c, ctx->jastrow.dtmp_c_device,
    sze * sizeof(double),
    0, 0,
    device_id, device_id
  );

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_dim_cord_vect(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.dim_cord_vect_date) {

    qmckl_exit_code rc =
      qmckl_compute_dim_cord_vect(context,
                                  ctx->jastrow.cord_num,
                                  &(ctx->jastrow.dim_cord_vect));
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.dim_cord_vect_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_cord_vect_full(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if dim_cord_vect is provided */
  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.cord_vect_full_date) {

    /* Allocate array */
    if (ctx->jastrow.cord_vect_full == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->jastrow.dim_cord_vect * ctx->nucleus.num * sizeof(double);
      double* cord_vect_full = (double*) qmckl_malloc(context, mem_info);

      if (cord_vect_full == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_cord_vect_full",
                               NULL);
      }
      ctx->jastrow.cord_vect_full = cord_vect_full;
    }

    rc = qmckl_compute_cord_vect_full(context,
                                      ctx->nucleus.num,
                                      ctx->jastrow.dim_cord_vect,
                                      ctx->jastrow.type_nucl_num,
                                      ctx->jastrow.type_nucl_vector,
                                      ctx->jastrow.cord_vector,
                                      ctx->jastrow.cord_vect_full);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.cord_vect_full_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_lkpm_combined_index(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if dim_cord_vect is provided */
  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.lkpm_combined_index_date) {

    /* Allocate array */
    if (ctx->jastrow.lkpm_combined_index == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->jastrow.dim_cord_vect * sizeof(int64_t);
      int64_t* lkpm_combined_index = (int64_t*) qmckl_malloc(context, mem_info);

      if (lkpm_combined_index == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_lkpm_combined_index",
                               NULL);
      }
      ctx->jastrow.lkpm_combined_index = lkpm_combined_index;
    }

    rc = qmckl_compute_lkpm_combined_index(context,
                                           ctx->jastrow.cord_num,
                                           ctx->jastrow.dim_cord_vect,
                                           ctx->jastrow.lkpm_combined_index);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.lkpm_combined_index_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_tmp_c(qmckl_context context)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if dim_cord_vect is provided */
  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.tmp_c_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.tmp_c);
      ctx->jastrow.tmp_c = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.tmp_c == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1)
                      * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);
      double* tmp_c = (double*) qmckl_malloc(context, mem_info);

      if (tmp_c == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_tmp_c",
                               NULL);
      }
      ctx->jastrow.tmp_c = tmp_c;
    }


    /* Choose the correct compute function (depending on offload type) */
#ifdef HAVE_HPC
    const bool gpu_offload = ctx->jastrow.gpu_offload;
#else
    const bool gpu_offload = false;
#endif

    if (gpu_offload) {
#ifdef HAVE_CUBLAS_OFFLOAD
      rc = qmckl_compute_tmp_c_cublas_offload(context,
                                              ctx->jastrow.cord_num,
                                              ctx->electron.num,
                                              ctx->nucleus.num,
                                              ctx->electron.walker.num,
                                              ctx->jastrow.een_rescaled_e,
                                              ctx->jastrow.een_rescaled_n,
                                              ctx->jastrow.tmp_c);
#elif HAVE_OPENACC_OFFLOAD
      rc = qmckl_compute_tmp_c_acc_offload(context,
                                           ctx->jastrow.cord_num,
                                           ctx->electron.num,
                                           ctx->nucleus.num,
                                           ctx->electron.walker.num,
                                           ctx->jastrow.een_rescaled_e,
                                           ctx->jastrow.een_rescaled_n,
                                           ctx->jastrow.tmp_c);
#elif HAVE_OPENMP_OFFLOAD
      rc = qmckl_compute_tmp_c_omp_offload(context,
                                           ctx->jastrow.cord_num,
                                           ctx->electron.num,
                                           ctx->nucleus.num,
                                           ctx->electron.walker.num,
                                           ctx->jastrow.een_rescaled_e,
                                           ctx->jastrow.een_rescaled_n,
                                           ctx->jastrow.tmp_c);
#else
      rc = QMCKL_FAILURE;
#endif
    } else {
      rc = qmckl_compute_tmp_c(context,
                               ctx->jastrow.cord_num,
                               ctx->electron.num,
                               ctx->nucleus.num,
                               ctx->electron.walker.num,
                               ctx->jastrow.een_rescaled_e,
                               ctx->jastrow.een_rescaled_n,
                               ctx->jastrow.tmp_c);
    }


    ctx->jastrow.tmp_c_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_dtmp_c(qmckl_context context)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if dim_cord_vect is provided */
  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.dtmp_c_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.dtmp_c);
      ctx->jastrow.dtmp_c = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.dtmp_c == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1)
                      * 4 * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);
      double* dtmp_c = (double*) qmckl_malloc(context, mem_info);

      if (dtmp_c == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_dtmp_c",
                               NULL);
      }
      ctx->jastrow.dtmp_c = dtmp_c;
    }


#ifdef HAVE_HPC
    const bool gpu_offload = ctx->jastrow.gpu_offload;
#else
    const bool gpu_offload = false;
#endif

    if (gpu_offload) {
#ifdef HAVE_CUBLAS_OFFLOAD
      rc = qmckl_compute_dtmp_c_cublas_offload(context,
                                            ctx->jastrow.cord_num,
                                            ctx->electron.num,
                                            ctx->nucleus.num,
                                            ctx->electron.walker.num,
                                            ctx->jastrow.een_rescaled_e_deriv_e,
                                            ctx->jastrow.een_rescaled_n,
                                            ctx->jastrow.dtmp_c);
#elif HAVE_OPENACC_OFFLOAD
      rc = qmckl_compute_dtmp_c_acc_offload(context,
                                            ctx->jastrow.cord_num,
                                            ctx->electron.num,
                                            ctx->nucleus.num,
                                            ctx->electron.walker.num,
                                            ctx->jastrow.een_rescaled_e_deriv_e,
                                            ctx->jastrow.een_rescaled_n,
                                            ctx->jastrow.dtmp_c);
#elif HAVE_OPENMP_OFFLOAD
      rc = qmckl_compute_dtmp_c_omp_offload(context,
                                            ctx->jastrow.cord_num,
                                            ctx->electron.num,
                                            ctx->nucleus.num,
                                            ctx->electron.walker.num,
                                            ctx->jastrow.een_rescaled_e_deriv_e,
                                            ctx->jastrow.een_rescaled_n,
                                            ctx->jastrow.dtmp_c);
#else
      rc = QMCKL_FAILURE;
#endif
    } else {
        rc = qmckl_compute_dtmp_c(context,
                                ctx->jastrow.cord_num,
                                ctx->electron.num,
                                ctx->nucleus.num,
                                ctx->electron.walker.num,
                                ctx->jastrow.een_rescaled_e_deriv_e,
                                ctx->jastrow.een_rescaled_n,
                                ctx->jastrow.dtmp_c);
    }

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }


    ctx->jastrow.dtmp_c_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_cord_vect_full_device(qmckl_context context, int device_id)
{


  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if dim_cord_vect is provided */

  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.cord_vect_full_date) {

    /* Allocate array */
    if (ctx->jastrow.cord_vect_full == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->jastrow.dim_cord_vect * ctx->nucleus.num * sizeof(double);
      double* cord_vect_full = (double*) qmckl_malloc_device(context, mem_info, device_id);

      if (cord_vect_full == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_cord_vect_full_device",
                               NULL);
      }
      ctx->jastrow.cord_vect_full_device = cord_vect_full;
    }

    rc = qmckl_compute_cord_vect_full_device(context,
                                             ctx->nucleus.num,
                                             ctx->jastrow.dim_cord_vect,
                                             ctx->jastrow.type_nucl_num,
                                             ctx->jastrow.type_nucl_vector_device,
                                             ctx->jastrow.cord_vector_device,
                                             ctx->jastrow.cord_vect_full_device,
                                             device_id);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.cord_vect_full_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_lkpm_combined_index_device(qmckl_context context, int device_id)
{


  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if dim_cord_vect is provided */
  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.lkpm_combined_index_date) {

    /* Allocate array */
    if (ctx->jastrow.lkpm_combined_index == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->jastrow.dim_cord_vect * sizeof(int64_t);
      int64_t* lkpm_combined_index = (int64_t*) qmckl_malloc_device(context, mem_info, device_id);

      if (lkpm_combined_index == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_lkpm_combined_index_device",
                               NULL);
      }
      ctx->jastrow.lkpm_combined_index_device = lkpm_combined_index;
    }

    rc = qmckl_compute_lkpm_combined_index_device(context,
                                                  ctx->jastrow.cord_num,
                                                  ctx->jastrow.dim_cord_vect,
                                                  ctx->jastrow.lkpm_combined_index_device,
                                                  device_id);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.lkpm_combined_index_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_tmp_c_device(qmckl_context context, int device_id)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if dim_cord_vect is provided */
  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.tmp_c_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      omp_target_free(ctx->jastrow.tmp_c_device, device_id);
      ctx->jastrow.tmp_c_device = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.tmp_c_device == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1)
                      * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);
      double* tmp_c = (double*) qmckl_malloc_device(context, mem_info, device_id);

      if (tmp_c == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_tmp_c_device",
                               NULL);
      }
      ctx->jastrow.tmp_c_device = tmp_c;
    }

      rc = qmckl_compute_tmp_c_device(context,
                                      ctx->jastrow.cord_num,
                                      ctx->electron.num,
                                      ctx->nucleus.num,
                                      ctx->electron.walker.num,
                                      ctx->jastrow.een_rescaled_e_device,
                                      ctx->jastrow.een_rescaled_n_device,
                                      ctx->jastrow.tmp_c_device,
                                      device_id);


    ctx->jastrow.tmp_c_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_dtmp_c_device(qmckl_context context, int device_id)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if dim_cord_vect is provided */
  qmckl_exit_code rc = qmckl_provide_dim_cord_vect(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.dtmp_c_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      omp_target_free(ctx->jastrow.dtmp_c_device, device_id);
      ctx->jastrow.dtmp_c_device = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.dtmp_c_device == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = (ctx->jastrow.cord_num) * (ctx->jastrow.cord_num + 1)
                      * 4 * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);
      double* dtmp_c = (double*) qmckl_malloc_device(context, mem_info, device_id);

      if (dtmp_c == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_dtmp_c_device",
                               NULL);
      }
      ctx->jastrow.dtmp_c_device = dtmp_c;
    }

    rc = qmckl_compute_dtmp_c_device(context,
                              ctx->jastrow.cord_num,
                              ctx->electron.num,
                              ctx->nucleus.num,
                              ctx->electron.walker.num,
                              ctx->jastrow.een_rescaled_e_deriv_e_device,
                              ctx->jastrow.een_rescaled_n_device,
                              ctx->jastrow.dtmp_c_device,
                              device_id);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }


    ctx->jastrow.dtmp_c_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_dim_cord_vect (
      const qmckl_context context,
      const int64_t cord_num,
      int64_t* const dim_cord_vect){

  int         lmax;


  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  *dim_cord_vect = 0;

  for (int p=2; p <= cord_num; ++p){
    for (int k=p-1; k >= 0; --k) {
      if (k != 0) {
        lmax = p - k;
      } else {
        lmax = p - k - 2;
      }
      for (int l = lmax; l >= 0; --l) {
        if ( ((p - k - l) & 1)==1) continue;
        *dim_cord_vect=*dim_cord_vect+1;
      }
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_cord_vect_full_hpc (
	  const qmckl_context context,
	  const int64_t nucl_num,
	  const int64_t dim_cord_vect,
	  const int64_t type_nucl_num,
	  const int64_t* type_nucl_vector,
	  const double* cord_vector,
	  double* const cord_vect_full ) {

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (nucl_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (type_nucl_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  if (dim_cord_vect <= 0) {
     return QMCKL_INVALID_ARG_5;
  }

  for (int i=0; i < dim_cord_vect; ++i) {
    for (int a=0; a < nucl_num; ++a){
      cord_vect_full[a + i*nucl_num] = cord_vector[(type_nucl_vector[a]-1)+i*type_nucl_num];
    }
  }

  return QMCKL_SUCCESS;
  }

qmckl_exit_code qmckl_compute_cord_vect_full (
	  const qmckl_context context,
	  const int64_t nucl_num,
	  const int64_t dim_cord_vect,
	  const int64_t type_nucl_num,
	  const int64_t* type_nucl_vector,
	  const double* cord_vector,
	  double* const cord_vect_full ) {

    #ifdef HAVE_HPC
      return qmckl_compute_cord_vect_full_hpc(context, nucl_num, dim_cord_vect, type_nucl_num, type_nucl_vector, cord_vector, cord_vect_full);
    #else
      return qmckl_compute_cord_vect_full_doc(context, nucl_num, dim_cord_vect, type_nucl_num, type_nucl_vector, cord_vector, cord_vect_full);
    #endif
    }

qmckl_exit_code qmckl_compute_cord_vect_full_device (
	  const qmckl_context context,
	  const int64_t nucl_num,
	  const int64_t dim_cord_vect,
	  const int64_t type_nucl_num,
	  const int64_t* type_nucl_vector,
	  const double* cord_vector,
	  double* const cord_vect_full,
      int device_id
 ) {

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (nucl_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (type_nucl_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  if (dim_cord_vect <= 0) {
     return QMCKL_INVALID_ARG_5;
  }

  #pragma omp target is_device_ptr(type_nucl_vector, cord_vector, cord_vect_full)
  {
  #pragma omp teams distribute parallel for collapse(2)
  for (int i=0; i < dim_cord_vect; ++i) {
    for (int a=0; a < nucl_num; ++a){
      cord_vect_full[a + i*nucl_num] = cord_vector[(type_nucl_vector[a]-1)+i*type_nucl_num];
    }
  }
  }
  return QMCKL_SUCCESS;
  }

qmckl_exit_code qmckl_compute_lkpm_combined_index (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_cord_vect,
      int64_t* const lkpm_combined_index ) {

  int kk, lmax, m;

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (dim_cord_vect <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

/*
*/
  kk = 0;
  for (int p = 2; p <= cord_num; ++p) {
    for (int k=(p-1); k >= 0; --k) {
      if (k != 0) {
        lmax = p - k;
      } else {
        lmax = p - k - 2;
      }
      for (int l=lmax; l >= 0; --l) {
        if (((p - k - l) & 1) == 1) continue;
        m = (p - k - l)/2;
        lkpm_combined_index[kk                  ] = l;
        lkpm_combined_index[kk +   dim_cord_vect] = k;
        lkpm_combined_index[kk + 2*dim_cord_vect] = p;
        lkpm_combined_index[kk + 3*dim_cord_vect] = m;
        kk = kk + 1;
      }
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_lkpm_combined_index_device (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_cord_vect,
      int64_t* const lkpm_combined_index,
      int device_id
) {

  int kk, lmax, m;

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (dim_cord_vect <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  #pragma omp target is_device_ptr(lkpm_combined_index)
  {
  kk = 0;
  #pragma omp parallel for collapse(2)
  for (int p = 2; p <= cord_num; ++p) {
    for (int k=(p-1); k >= 0; --k) {
      if (k != 0) {
        lmax = p - k;
      } else {
        lmax = p - k - 2;
      }
      for (int l=lmax; l >= 0; --l) {
        if (((p - k - l) & 1) == 1) continue;
        m = (p - k - l)/2;
        lkpm_combined_index[kk                  ] = l;
        lkpm_combined_index[kk +   dim_cord_vect] = k;
        lkpm_combined_index[kk + 2*dim_cord_vect] = p;
        lkpm_combined_index[kk + 3*dim_cord_vect] = m;
        kk = kk + 1;
      }
    }
  }
  }

  return QMCKL_SUCCESS;
}

/* Doc */

/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_tmp_c */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*     #+NAME: qmckl_factor_tmp_c_args */
/*     | Variable         | Type                                                             | In/Out | Description                       | */
/*     |------------------+------------------------------------------------------------------+--------+-----------------------------------| */
/*     | ~context~        | ~qmckl_context~                                                  | in     | Global state                      | */
/*     | ~cord_num~       | ~int64_t~                                                        | in     | Order of polynomials              | */
/*     | ~elec_num~       | ~int64_t~                                                        | in     | Number of electrons               | */
/*     | ~nucl_num~       | ~int64_t~                                                        | in     | Number of nucleii                 | */
/*     | ~walk_num~       | ~int64_t~                                                        | in     | Number of walkers                 | */
/*     | ~een_rescaled_e~ | ~double[walk_num][0:cord_num][elec_num][elec_num]~               | in     | Electron-electron rescaled factor | */
/*     | ~een_rescaled_n~ | ~double[walk_num][0:cord_num][nucl_num][elec_num]~               | in     | Electron-nucleus rescaled factor  | */
/*     | ~tmp_c~          | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~ | out    | vector of non-zero coefficients   | */


qmckl_exit_code qmckl_compute_tmp_c (const qmckl_context context,
                                     const int64_t cord_num,
                                     const int64_t elec_num,
                                     const int64_t nucl_num,
                                     const int64_t walk_num,
                                     const double* een_rescaled_e,
                                     const double* een_rescaled_n,
                                     double* const tmp_c )
{
#ifdef HAVE_HPC
  return qmckl_compute_tmp_c_hpc(context, cord_num, elec_num, nucl_num, walk_num, een_rescaled_e, een_rescaled_n, tmp_c);
#else
  return qmckl_compute_tmp_c_doc(context, cord_num, elec_num, nucl_num, walk_num, een_rescaled_e, een_rescaled_n, tmp_c);
#endif
}

/* CPU                                                           :noexport: */


qmckl_exit_code qmckl_compute_tmp_c_hpc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e,
      const double* een_rescaled_n,
      double* const tmp_c ) {

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  if (walk_num <= 0) {
    return QMCKL_INVALID_ARG_5;
  }

  qmckl_exit_code info = QMCKL_SUCCESS;

  const char  TransA = 'N';
  const char  TransB = 'N';
  const double alpha = 1.0;
  const double beta  = 0.0;

  const int64_t M = elec_num;
  const int64_t N = nucl_num*(cord_num + 1);
  const int64_t K = elec_num;

  const int64_t LDA = elec_num;
  const int64_t LDB = elec_num;
  const int64_t LDC = elec_num;

  const int64_t af = elec_num*elec_num;
  const int64_t bf = elec_num*nucl_num*(cord_num+1);
  const int64_t cf = bf;

#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int64_t nw=0; nw < walk_num; ++nw) {
    for (int64_t i=0; i<cord_num; ++i){
      info = qmckl_dgemm(context, TransA, TransB, M, N, K, alpha,
                         &(een_rescaled_e[af*(i+nw*(cord_num+1))]), LDA,
                         &(een_rescaled_n[bf*nw]), LDB, beta,
                         &(tmp_c[cf*(i+nw*cord_num)]), LDC);
    }
  }

  return info;
}

/* OpenACC offload                                               :noexport: */


#ifdef HAVE_OPENACC_OFFLOAD
qmckl_exit_code
qmckl_compute_tmp_c_acc_offload (const qmckl_context context,
                                 const int64_t cord_num,
                                 const int64_t elec_num,
                                 const int64_t nucl_num,
                                 const int64_t walk_num,
                                 const double* een_rescaled_e,
                                 const double* een_rescaled_n,
                                 double* const tmp_c )
{

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  // Compute array access strides:
  // For tmp_c...
  const int64_t stride_k_c  = elec_num;
  const int64_t stride_j_c  = stride_k_c * nucl_num;
  const int64_t stride_i_c  = stride_j_c * (cord_num+1);
  const int64_t stride_nw_c = stride_i_c * cord_num;
  // For een_rescaled_e...
  const int64_t stride_m_e  = elec_num;
  const int64_t stride_i_e  = stride_m_e * elec_num;
  const int64_t stride_nw_e = stride_i_e * (cord_num+1);
  // For een_rescaled_n...
  const int64_t stride_k_n  = elec_num;
  const int64_t stride_j_n  = stride_k_n * nucl_num;
  const int64_t stride_nw_n = stride_j_n * (cord_num+1);

  const int64_t size_tmp_c = elec_num*nucl_num*(cord_num+1)*cord_num*walk_num;
  const int64_t size_e = walk_num*(cord_num+1)*elec_num*elec_num;
  const int64_t size_n = walk_num*(cord_num+1)*nucl_num*elec_num;

  #pragma acc parallel copyout(tmp_c [0:size_tmp_c]) copyin(een_rescaled_e[0:size_e], een_rescaled_n[0:size_n])
  {
  #pragma acc loop independent gang worker vector collapse(4)
  for (int nw=0; nw < walk_num; ++nw) {
    for (int i=0; i<cord_num; ++i){

      // Replacement for single DGEMM
      for (int jk=0; jk<nucl_num*(cord_num+1); jk++) {
        for (int l=0; l<elec_num; l++) {

		  int index_e_base = l + i*stride_i_e + nw*stride_nw_e;
		  int index_n_base = jk*stride_k_n + nw*stride_nw_n;

          // Single reduction
		  double sum = 0.;
          for (int m=0; m<elec_num; m++) {
            sum +=
            een_rescaled_e[index_e_base + m*stride_m_e] *
            een_rescaled_n[index_n_base + m];
          }
          tmp_c[l + jk*stride_k_c + i*stride_i_c + nw*stride_nw_c] = sum;

        }
      }
    }
  }
  }

  return QMCKL_SUCCESS;
}
#endif

/* OpenMP offload                                                :noexport: */


#ifdef HAVE_OPENMP_OFFLOAD
qmckl_exit_code
qmckl_compute_tmp_c_omp_offload (const qmckl_context context,
                                 const int64_t cord_num,
                                 const int64_t elec_num,
                                 const int64_t nucl_num,
                                 const int64_t walk_num,
                                 const double* een_rescaled_e,
                                 const double* een_rescaled_n,
                                 double* const tmp_c )
{

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  // Compute array access strides:
  // For tmp_c...
  const int64_t stride_k_c  = elec_num;
  const int64_t stride_j_c  = stride_k_c * nucl_num;
  const int64_t stride_i_c  = stride_j_c * (cord_num+1);
  const int64_t stride_nw_c = stride_i_c * cord_num;
  // For een_rescaled_e...
  const int64_t stride_m_e  = elec_num;
  const int64_t stride_i_e  = stride_m_e * elec_num;
  const int64_t stride_nw_e = stride_i_e * (cord_num+1);
  // For een_rescaled_n...
  const int64_t stride_k_n  = elec_num;
  const int64_t stride_j_n  = stride_k_n * nucl_num;
  const int64_t stride_nw_n = stride_j_n * (cord_num+1);

  const int64_t size_tmp_c = elec_num*nucl_num*(cord_num+1)*cord_num*walk_num;
  const int64_t size_e = walk_num*(cord_num+1)*elec_num*elec_num;
  const int64_t size_n = walk_num*(cord_num+1)*nucl_num*elec_num;

  #pragma omp target data map(from:tmp_c[0:size_tmp_c]) map(to:een_rescaled_e[0:size_e], een_rescaled_n[0:size_n])
  {
  #pragma omp target teams distribute parallel for collapse(4)
  for (int nw=0; nw < walk_num; ++nw) {
    for (int i=0; i<cord_num; ++i){

      // Replacement for single DGEMM
      for (int jk=0; jk<nucl_num*(cord_num+1); jk++) {
        for (int l=0; l<elec_num; l++) {

		  int index_e_base = l + i*stride_i_e + nw*stride_nw_e;
		  int index_n_base = jk*stride_k_n + nw*stride_nw_n;

          // Single reduction
		  double sum = 0.;
		  #pragma omp simd reduction(+:sum)
          for (int m=0; m<elec_num; m++) {
            sum +=
            een_rescaled_e[index_e_base + m*stride_m_e] *
            een_rescaled_n[index_n_base + m];
          }
          tmp_c[l + jk*stride_k_c + i*stride_i_c + nw*stride_nw_c] = sum;

        }
      }
    }
  }
  }

  return QMCKL_SUCCESS;
}
#endif

/* cuBLAS offload                                                :noexport: */


#ifdef HAVE_CUBLAS_OFFLOAD
qmckl_exit_code
qmckl_compute_tmp_c_cublas_offload (const qmckl_context context,
                                    const int64_t cord_num,
                                    const int64_t elec_num,
                                    const int64_t nucl_num,
                                    const int64_t walk_num,
                                    const double* een_rescaled_e,
                                    const double* een_rescaled_n,
                                    double* const tmp_c )
{
  qmckl_exit_code info;

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  //cuBLAS initialization
  cublasHandle_t handle;
  if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS)
  {
    fprintf(stdout, "CUBLAS initialization failed!\n");
    exit(EXIT_FAILURE);
  }

  const double alpha = 1.0;
  const double beta  = 0.0;

  const int64_t M = elec_num;
  const int64_t N = nucl_num*(cord_num + 1);
  const int64_t K = elec_num;

  const int64_t LDA = elec_num;
  const int64_t LDB = elec_num;
  const int64_t LDC = elec_num;

  const int64_t af = elec_num*elec_num;
  const int64_t bf = elec_num*nucl_num*(cord_num+1);
  const int64_t cf = bf;

  #pragma omp target enter data map(to:een_rescaled_e[0:elec_num*elec_num*(cord_num+1)*walk_num],een_rescaled_n[0:M*N*walk_num],tmp_c[0:elec_num*nucl_num*(cord_num+1)*cord_num*walk_num])
  #pragma omp target data use_device_ptr(een_rescaled_e,een_rescaled_n,tmp_c)
  {
  for (int nw=0; nw < walk_num; ++nw) {

    int cublasError = cublasDgemmStridedBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, K, &alpha,
                                    &(een_rescaled_e[nw*(cord_num+1)]),
                                    LDA, af,
                                    &(een_rescaled_n[bf*nw]),
                                    LDB, 0,
                                    &beta,
                                    &(tmp_c[nw*cord_num]),
                                    LDC, cf, cord_num);
  }
  }
  #pragma omp target exit data map(from:tmp_c[0:elec_num*nucl_num*(cord_num+1)*cord_num*walk_num])

  cublasDestroy(handle);
  return info;
  }
#endif

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_compute_tmp_c_device (const qmckl_context context,
                            const int64_t cord_num,
                            const int64_t elec_num,
                            const int64_t nucl_num,
                            const int64_t walk_num,
                            const double* een_rescaled_e,
                            const double* een_rescaled_n,
                            double* const tmp_c,

                            int device_id )
{

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  // Compute array access strides:
  // For tmp_c...
  const int64_t stride_k_c  = elec_num;
  const int64_t stride_j_c  = stride_k_c * nucl_num;
  const int64_t stride_i_c  = stride_j_c * (cord_num+1);
  const int64_t stride_nw_c = stride_i_c * cord_num;
  // For een_rescaled_e...
  const int64_t stride_m_e  = elec_num;
  const int64_t stride_i_e  = stride_m_e * elec_num;
  const int64_t stride_nw_e = stride_i_e * (cord_num+1);
  // For een_rescaled_n...
  const int64_t stride_k_n  = elec_num;
  const int64_t stride_j_n  = stride_k_n * nucl_num;
  const int64_t stride_nw_n = stride_j_n * (cord_num+1);

  const int64_t size_tmp_c = elec_num*nucl_num*(cord_num+1)*cord_num*walk_num;
  const int64_t size_e = walk_num*(cord_num+1)*elec_num*elec_num;
  const int64_t size_n = walk_num*(cord_num+1)*nucl_num*elec_num;

  #pragma omp target is_device_ptr(tmp_c, een_rescaled_e, een_rescaled_n)
  {
  #pragma omp teams distribute parallel for collapse(4)
  for (int nw=0; nw < walk_num; ++nw) {
    for (int i=0; i<cord_num; ++i){

      // Replacement for single DGEMM
      for (int jk=0; jk<nucl_num*(cord_num+1); jk++) {
        for (int l=0; l<elec_num; l++) {

          int index_e_base = l + i*stride_i_e + nw*stride_nw_e;
          int index_n_base = jk*stride_k_n + nw*stride_nw_n;

          // Single reduction
          double sum = 0.;
          #pragma omp simd reduction(+:sum)
          for (int m=0; m<elec_num; m++) {
            sum +=
            een_rescaled_e[index_e_base + m*stride_m_e] *
            een_rescaled_n[index_n_base + m];
          }
          tmp_c[l + jk*stride_k_c + i*stride_i_c + nw*stride_nw_c] = sum;

        }
      }
    }
  }
  }

  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code
qmckl_compute_dtmp_c (const qmckl_context context,
                      const int64_t cord_num,
                      const int64_t elec_num,
                      const int64_t nucl_num,
                      const int64_t walk_num,
                      const double* een_rescaled_e_deriv_e,
                      const double* een_rescaled_n,
                      double* const dtmp_c )
{
#ifdef HAVE_HPC
  return qmckl_compute_dtmp_c_hpc (context, cord_num, elec_num, nucl_num, walk_num, een_rescaled_e_deriv_e,
                                   een_rescaled_n, dtmp_c );
#else
  return qmckl_compute_dtmp_c_doc (context, cord_num, elec_num, nucl_num, walk_num, een_rescaled_e_deriv_e,
                                   een_rescaled_n, dtmp_c );
#endif
}

/* CPU                                                           :noexport: */


qmckl_exit_code
qmckl_compute_dtmp_c_hpc (const qmckl_context context,
                          const int64_t cord_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t walk_num,
                          const double* een_rescaled_e_deriv_e,
                          const double* een_rescaled_n,
                          double* const dtmp_c )
{

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  if (walk_num <= 0) {
     return QMCKL_INVALID_ARG_5;
  }

  qmckl_exit_code  info = QMCKL_SUCCESS;

  const char  TransA = 'N';
  const char  TransB = 'N';
  const double alpha = 1.0;
  const double beta  = 0.0;

  const int64_t M = 4*elec_num;
  const int64_t N = nucl_num*(cord_num + 1);
  const int64_t K = elec_num;

  const int64_t LDA = 4*elec_num;
  const int64_t LDB = elec_num;
  const int64_t LDC = 4*elec_num;

  const int64_t af = elec_num*elec_num*4;
  const int64_t bf = elec_num*nucl_num*(cord_num+1);
  const int64_t cf = elec_num*4*nucl_num*(cord_num+1);

#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int64_t nw=0; nw < walk_num; ++nw) {
    for (int64_t i=0; i < cord_num; ++i) {
      info = qmckl_dgemm(context, TransA, TransB, M, N, K, alpha,
                         &(een_rescaled_e_deriv_e[af*(i+nw*(cord_num+1))]), LDA,
                         &(een_rescaled_n[bf*nw]), LDB, beta,
                         &(dtmp_c[cf*(i+nw*cord_num)]), LDC);
    }
  }

  return info;
}

/* OpenACC offload                                               :noexport: */


#ifdef HAVE_OPENACC_OFFLOAD
qmckl_exit_code
qmckl_compute_dtmp_c_acc_offload (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c ) {

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  // Compute strides...
  // For dtmp_c
  const int64_t stride_l_d  = elec_num;
  const int64_t stride_k_d  = stride_l_d * 4;
  const int64_t stride_j_d  = stride_k_d * nucl_num;
  const int64_t stride_i_d  = stride_j_d * (cord_num+1);
  const int64_t stride_nw_d = stride_i_d * cord_num;
  // For een_rescaled_e_deriv_e
  const int64_t stride_l_e  = elec_num;
  const int64_t stride_n_e  = stride_l_e * 4;
  const int64_t stride_i_e  = stride_n_e * elec_num;
  const int64_t stride_nw_e = stride_i_e * cord_num;
  // For een_rescaled_n
  const int64_t stride_k_n  = elec_num;
  const int64_t stride_j_n  = stride_k_n * nucl_num;
  const int64_t stride_nw_n = stride_j_n * (cord_num+1);

  const int64_t size_dtmp_c = walk_num*cord_num*(cord_num+1)*nucl_num*4*elec_num;
  const int64_t size_n = walk_num*(cord_num+1)*nucl_num*elec_num;
  const int64_t size_e = walk_num*(cord_num+1)*elec_num*4*elec_num;

  #pragma acc parallel copyout(dtmp_c [0:size_dtmp_c]) copyin(een_rescaled_e_deriv_e[0:size_e], een_rescaled_n[0:size_n])
  {
  #pragma acc loop independent gang worker vector collapse(4)
  for (int nw=0; nw < walk_num; nw++) {
    for (int i=0; i < cord_num; i++) {

      // Single DGEMM
      for(int jk=0; jk<nucl_num*(cord_num+1); jk++) {
        for(int ml=0; ml<4*elec_num; ml++) {

          // Single reduction
		  int index_n_base = jk * stride_k_n  + nw * stride_nw_n;
		  int index_e_base = ml + i * stride_i_e + nw * stride_nw_e;
		  double sum = 0.;
          for(int n=0; n<elec_num; n++){
            sum +=
            een_rescaled_e_deriv_e[index_e_base + n * stride_n_e] *
            een_rescaled_n[index_n_base + n];
          }
          dtmp_c[ml + jk * stride_k_d + i * stride_i_d + nw * stride_nw_d] = sum;

        }
      }
    }
  }
  }

  return QMCKL_SUCCESS;
}

#endif

/* OpenMP offload                                                :noexport: */


#ifdef HAVE_OPENMP_OFFLOAD
qmckl_exit_code qmckl_compute_dtmp_c_omp_offload (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c ) {

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  // Compute strides...
  // For dtmp_c
  const int64_t stride_l_d  = elec_num;
  const int64_t stride_k_d  = stride_l_d * 4;
  const int64_t stride_j_d  = stride_k_d * nucl_num;
  const int64_t stride_i_d  = stride_j_d * (cord_num+1);
  const int64_t stride_nw_d = stride_i_d * cord_num;
  // For een_rescaled_e_deriv_e
  const int64_t stride_l_e  = elec_num;
  const int64_t stride_n_e  = stride_l_e * 4;
  const int64_t stride_i_e  = stride_n_e * elec_num;
  const int64_t stride_nw_e = stride_i_e * cord_num;
  // For een_rescaled_n
  const int64_t stride_k_n  = elec_num;
  const int64_t stride_j_n  = stride_k_n * nucl_num;
  const int64_t stride_nw_n = stride_j_n * (cord_num+1);

  const int64_t size_dtmp_c = walk_num*cord_num*(cord_num+1)*nucl_num*4*elec_num;
  const int64_t size_n = walk_num*(cord_num+1)*nucl_num*elec_num;
  const int64_t size_e = walk_num*(cord_num+1)*elec_num*4*elec_num;


  double sum = 0.;
  #pragma omp target data map(from:dtmp_c[0:size_dtmp_c]) map(to:een_rescaled_e_deriv_e[0:size_e], een_rescaled_n[0:size_n])
  {
  #pragma omp target teams distribute parallel for collapse(4)
  for (int nw=0; nw < walk_num; nw++) {
    for (int i=0; i < cord_num; i++) {

      // Single DGEMM
      for(int jk=0; jk<nucl_num*(cord_num+1); jk++) {
        for(int ml=0; ml<4*elec_num; ml++) {

          // Single reduction
          int index_n_base = jk * stride_k_n  + nw * stride_nw_n;
          int index_e_base = ml + i * stride_i_e + nw * stride_nw_e;
          sum = 0.;
          #pragma omp simd reduction(+:sum)
          for(int n=0; n<elec_num; n++){
            sum +=
            een_rescaled_e_deriv_e[index_e_base + n * stride_n_e] *
            een_rescaled_n[index_n_base + n];
          }
          dtmp_c[ml + jk * stride_k_d + i * stride_i_d + nw * stride_nw_d] = sum;

        }
      }
    }
  }
  }

  return QMCKL_SUCCESS;
}
#endif

/* cuBLAS offload                                                :noexport: */


#ifdef HAVE_CUBLAS_OFFLOAD
qmckl_exit_code
qmckl_compute_dtmp_c_cublas_offload (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c ) {

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  if (walk_num <= 0) {
    return QMCKL_INVALID_ARG_5;
  }

  qmckl_exit_code  info = QMCKL_SUCCESS;

  //cuBLAS initialization
  cublasHandle_t handle;
  if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS)
  {
    fprintf(stdout, "CUBLAS initialization failed!\n");
    exit(EXIT_FAILURE);
  }

  const double alpha = 1.0;
  const double beta  = 0.0;

  const int64_t M = 4*elec_num;
  const int64_t N = nucl_num*(cord_num + 1);
  const int64_t K = elec_num;

  const int64_t LDA = 4*elec_num;
  const int64_t LDB = elec_num;
  const int64_t LDC = 4*elec_num;

  const int64_t af = elec_num*elec_num*4;
  const int64_t bf = elec_num*nucl_num*(cord_num+1);
  const int64_t cf = elec_num*4*nucl_num*(cord_num+1);

  #pragma omp target enter data map(to:een_rescaled_e_deriv_e[0:elec_num*4*elec_num*(cord_num+1)*walk_num], een_rescaled_n[0:elec_num*nucl_num*(cord_num+1)*walk_num], dtmp_c[0:elec_num*4*nucl_num*(cord_num+1)*cord_num*walk_num])
  #pragma omp target data use_device_ptr(een_rescaled_e_deriv_e, een_rescaled_n, dtmp_c)
  {
  #pragma omp parallel for
  for (int64_t nw=0; nw < walk_num; ++nw) {
    int cublasError = cublasDgemmStridedBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, K, &alpha,
                                      &(een_rescaled_e_deriv_e[(nw*(cord_num+1))]),
                                      LDA, af,
                                      &(een_rescaled_n[bf*nw]), LDB, 0,
                                      &beta,
                                      &(dtmp_c[(nw*cord_num)]),
                                      LDC, cf, cord_num);

  }
  }

  #pragma omp target exit data map(from:dtmp_c[0:cf*cord_num*walk_num])

  cublasDestroy(handle);
  return info;
}
#endif

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_compute_dtmp_c_device (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_deriv_e,
      const double* een_rescaled_n,
      double* const dtmp_c ) {

  if (context == QMCKL_NULL_CONTEXT) {
     return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num <= 0) {
     return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
     return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
     return QMCKL_INVALID_ARG_4;
  }

  // Compute strides...
  // For dtmp_c
  const int64_t stride_l_d  = elec_num;
  const int64_t stride_k_d  = stride_l_d * 4;
  const int64_t stride_j_d  = stride_k_d * nucl_num;
  const int64_t stride_i_d  = stride_j_d * (cord_num+1);
  const int64_t stride_nw_d = stride_i_d * cord_num;
  // For een_rescaled_e_deriv_e
  const int64_t stride_l_e  = elec_num;
  const int64_t stride_n_e  = stride_l_e * 4;
  const int64_t stride_i_e  = stride_n_e * elec_num;
  const int64_t stride_nw_e = stride_i_e * cord_num;
  // For een_rescaled_n
  const int64_t stride_k_n  = elec_num;
  const int64_t stride_j_n  = stride_k_n * nucl_num;
  const int64_t stride_nw_n = stride_j_n * (cord_num+1);

  const int64_t size_dtmp_c = walk_num*cord_num*(cord_num+1)*nucl_num*4*elec_num;
  const int64_t size_n = walk_num*(cord_num+1)*nucl_num*elec_num;
  const int64_t size_e = walk_num*(cord_num+1)*elec_num*4*elec_num;


  double sum = 0.;
  #pragma omp target is_device_ptr(dtmp_c, een_rescaled_e_deriv_e, een_rescaled_n)
  {
  #pragma omp teams distribute parallel for collapse(4)
  for (int nw=0; nw < walk_num; nw++) {
    for (int i=0; i < cord_num; i++) {

      // Single DGEMM
      for(int jk=0; jk<nucl_num*(cord_num+1); jk++) {
        for(int ml=0; ml<4*elec_num; ml++) {

          // Single reduction
          int index_n_base = jk * stride_k_n  + nw * stride_nw_n;
          int index_e_base = ml + i * stride_i_e + nw * stride_nw_e;
          sum = 0.;
          #pragma omp simd reduction(+:sum)
          for(int n=0; n<elec_num; n++){
            sum +=
            een_rescaled_e_deriv_e[index_e_base + n * stride_n_e] *
            een_rescaled_n[index_n_base + n];
          }
          dtmp_c[ml + jk * stride_k_d + i * stride_i_d + nw * stride_nw_d] = sum;

        }
      }
    }
  }
  }

  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code
qmckl_get_jastrow_factor_een(qmckl_context context,
                             double* const factor_een,
                             const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_factor_een(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_een",
                           "Array too small. Expected walk_num");
  }
  memcpy(factor_een, ctx->jastrow.factor_een, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_factor_een(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_een_rescaled_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_een_rescaled_n(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_cord_vect_full(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_lkpm_combined_index(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if tmp_c is provided */
  rc = qmckl_provide_tmp_c(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.factor_een_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.factor_een);
      ctx->jastrow.factor_een = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.factor_een == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* factor_een = (double*) qmckl_malloc(context, mem_info);

      if (factor_een == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_factor_een",
                               NULL);
      }
      ctx->jastrow.factor_een = factor_een;
    }

    rc = qmckl_compute_factor_een(context,
                                  ctx->electron.walker.num,
                                  ctx->electron.num,
                                  ctx->nucleus.num,
                                  ctx->jastrow.cord_num,
                                  ctx->jastrow.dim_cord_vect,
                                  ctx->jastrow.cord_vect_full,
                                  ctx->jastrow.lkpm_combined_index,
                                  ctx->jastrow.tmp_c,
                                  ctx->jastrow.een_rescaled_n,
                                  ctx->jastrow.factor_een);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.factor_een_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_factor_een_deriv_e(qmckl_context context,
                                     double* const factor_een_deriv_e,
                                     const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_factor_een_deriv_e(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_een_deriv_e",
                           "Array too small. Expected 4*walk_num*elec_num");
  }
  memcpy(factor_een_deriv_e, ctx->jastrow.factor_een_deriv_e, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_get_jastrow_factor_een_deriv_e_device(
  qmckl_context context,
  double* const factor_een_deriv_e,
  const int64_t size_max,
  int device_id
)
{


  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_factor_een_deriv_e_device(context, device_id);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_factor_een_deriv_e_device",
                           "Array too small. Expected 4*walk_num*elec_num");
  }
  // memcpy(factor_een_deriv_e, ctx->jastrow.factor_een_deriv_e, sze*sizeof(double));
  omp_target_memcpy(factor_een_deriv_e, ctx->jastrow.factor_een_deriv_e_device,
                    sze*sizeof(double),
                    0, 0,
                    device_id, device_id
                    );


  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code qmckl_provide_factor_een_deriv_e(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_een_rescaled_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_een_rescaled_n(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_een_rescaled_e_deriv_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_een_rescaled_n_deriv_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_cord_vect_full(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_lkpm_combined_index(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if tmp_c is provided */
  rc = qmckl_provide_tmp_c(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if dtmp_c is provided */
  rc = qmckl_provide_dtmp_c(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.factor_een_deriv_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->jastrow.factor_een_deriv_e);
      ctx->jastrow.factor_een_deriv_e = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.factor_een_deriv_e == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * ctx->electron.walker.num * sizeof(double);
      double* factor_een_deriv_e = (double*) qmckl_malloc(context, mem_info);

      if (factor_een_deriv_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_factor_een_deriv_e",
                               NULL);
      }
      ctx->jastrow.factor_een_deriv_e = factor_een_deriv_e;
    }

    rc = qmckl_compute_factor_een_deriv_e(context,
                                          ctx->electron.walker.num,
                                          ctx->electron.num,
                                          ctx->nucleus.num,
                                          ctx->jastrow.cord_num,
                                          ctx->jastrow.dim_cord_vect,
                                          ctx->jastrow.cord_vect_full,
                                          ctx->jastrow.lkpm_combined_index,
                                          ctx->jastrow.tmp_c,
                                          ctx->jastrow.dtmp_c,
                                          ctx->jastrow.een_rescaled_n,
                                          ctx->jastrow.een_rescaled_n_deriv_e,
                                          ctx->jastrow.factor_een_deriv_e);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.factor_een_deriv_e_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_provide_factor_een_deriv_e_device(qmckl_context context, int device_id)
{
  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_een_rescaled_e_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_een_rescaled_n_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_een_rescaled_e_deriv_e_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_een_rescaled_n_deriv_e_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_cord_vect_full_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_lkpm_combined_index_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if tmp_c is provided */
  rc = qmckl_provide_tmp_c_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if dtmp_c is provided */
  rc = qmckl_provide_dtmp_c_device(context, device_id);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow.factor_een_deriv_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      omp_target_free(ctx->jastrow.factor_een_deriv_e_device, device_id);
      ctx->jastrow.factor_een_deriv_e_device = NULL;
    }

    /* Allocate array */
    if (ctx->jastrow.factor_een_deriv_e_device == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * ctx->electron.walker.num * sizeof(double);
      double* factor_een_deriv_e = (double*) qmckl_malloc_device(context, mem_info, device_id);

      if (factor_een_deriv_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_factor_een_deriv_e_device",
                               NULL);
      }
      ctx->jastrow.factor_een_deriv_e_device = factor_een_deriv_e;
    }

    rc = qmckl_compute_factor_een_deriv_e_device(context,
                                          ctx->electron.walker.num,
                                          ctx->electron.num,
                                          ctx->nucleus.num,
                                          ctx->jastrow.cord_num,
                                          ctx->jastrow.dim_cord_vect,
                                          ctx->jastrow.cord_vect_full_device,
                                          ctx->jastrow.lkpm_combined_index_device,
                                          ctx->jastrow.tmp_c_device,
                                          ctx->jastrow.dtmp_c_device,
                                          ctx->jastrow.een_rescaled_n_device,
                                          ctx->jastrow.een_rescaled_n_deriv_e_device,
                                          ctx->jastrow.factor_een_deriv_e_device,
                                          device_id);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow.factor_een_deriv_e_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}
#endif

#ifdef HAVE_DEVICE_POINTERS
    qmckl_exit_code qmckl_compute_factor_een_deriv_e_device (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const int64_t dim_cord_vect,
      const double* cord_vect_full,
      const int64_t* lkpm_combined_index,
      const double* tmp_c,
      const double* dtmp_c,
      const double* een_rescaled_n,
      const double* een_rescaled_n_deriv_e,
      double* const factor_een_deriv_e,
      int device_id ){


      qmckl_exit_code rc = QMCKL_SUCCESS;

      int l,k,p,m;
      double accu, accu2, cn;

      if(context == QMCKL_NULL_CONTEXT){
        return QMCKL_INVALID_CONTEXT;
      }

      if(walk_num < 0){
        return QMCKL_INVALID_ARG_2;
      }

      if(elec_num < 0){
        return QMCKL_INVALID_ARG_3;
      }

      if(nucl_num < 0){
        return QMCKL_INVALID_ARG_4;
      }

      if(cord_num < 0){
        return QMCKL_INVALID_ARG_5;
      }

/*
  integer*8             , intent(in)  :: lkpm_combined_index(dim_cord_vect,4)
  double precision      , intent(in)  :: cord_vect_full(nucl_num, dim_cord_vect)
  double precision      , intent(in)  :: tmp_c(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  double precision      , intent(in)  :: dtmp_c(elec_num, 4, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  double precision      , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  double precision      , intent(in)  :: een_rescaled_n_deriv_e(elec_num, 4, nucl_num, 0:cord_num, walk_num)
  double precision      , intent(out) :: factor_een_deriv_e(elec_num,4,walk_num)
*/

      // Pre-compute arrays offsets
      // For tmp_c:
      int tmp_offset_1 = elec_num;
      int tmp_offset_2 = tmp_offset_1 * nucl_num;
      int tmp_offset_3 = tmp_offset_1 * (cord_num+1);
      int tmp_offset_4 = tmp_offset_1 * cord_num;

      // For dtmp_c:
      int dtmp_offset_1 = elec_num;
      int dtmp_offset_2 = dtmp_offset_1 * 4;
      int dtmp_offset_3 = dtmp_offset_2 * nucl_num;
      int dtmp_offset_4 = dtmp_offset_3 * (cord_num+1);
      int dtmp_offset_5 = dtmp_offset_4 * cord_num;

      // For een_rescaled_n:
      int een_r_offset_1 = elec_num;
      int een_r_offset_2 = een_r_offset_1 * nucl_num;
      int een_r_offset_3 = een_r_offset_2 * (cord_num+1);

      // For een_rescaled_n_deriv_e:
      int een_rd_offset_1 = elec_num;
      int een_rd_offset_2 = een_rd_offset_1 * 4;
      int een_rd_offset_3 = een_rd_offset_2 * nucl_num;
      int een_rd_offset_4 = een_rd_offset_3 * (cord_num+1);

      // For factor_een_deriv_e:
      int factor_offset_1 = elec_num;
      int factor_offset_2 = factor_offset_1 * 4;


      #pragma omp target is_device_ptr(cord_vect_full, lkpm_combined_index, tmp_c, dtmp_c, een_rescaled_n, een_rescaled_n_deriv_e, factor_een_deriv_e)
      {

      #pragma omp parallel for collapse(2)
      for(int nw = 0; nw < walk_num; nw++){
        for(int n =0; n < dim_cord_vect; n++){


          l = lkpm_combined_index[n + 1*dim_cord_vect];
          k = lkpm_combined_index[n + 2*dim_cord_vect];
          p = lkpm_combined_index[n + 3*dim_cord_vect];
          m = lkpm_combined_index[n + 4*dim_cord_vect];

          for(int a = 1; a<nucl_num; a++) {
            cn = cord_vect_full[a + n*nucl_num];
            if(cn == 0.0) continue;

            for(int ii = 1; ii<4; ii++) {
              for(int j = 1; j<elec_num; j++) {
                factor_een_deriv_e[j + ii*factor_offset_1 + nw*factor_offset_2] = factor_een_deriv_e[j + ii*factor_offset_1 + nw*factor_offset_2] + (tmp_c[j + a*tmp_offset_1 + m*tmp_offset_2 + k*tmp_offset_3 + nw*tmp_offset_4] * een_rescaled_n_deriv_e[j + ii*een_rd_offset_1 + a*een_rd_offset_2 + (m+l)*een_rd_offset_3 + nw*een_rd_offset_4] + dtmp_c[j + ii*dtmp_offset_1 + a*dtmp_offset_2 + m*dtmp_offset_3 + k*dtmp_offset_4 + nw*dtmp_offset_5] * een_rescaled_n[j + a*een_r_offset_1 + (m+l)*een_r_offset_2 + nw*een_r_offset_3] + dtmp_c[j + ii*dtmp_offset_1 + a*dtmp_offset_2 + (m+l)*dtmp_offset_3 + k*dtmp_offset_4 + nw*dtmp_offset_5] * een_rescaled_n[j + a*een_r_offset_1 + m*een_r_offset_2 + nw*een_r_offset_3] + tmp_c[j + a*tmp_offset_1 + (m+l)*tmp_offset_2 + k*tmp_offset_3 + nw*tmp_offset_4] * een_rescaled_n_deriv_e[j + ii*een_rd_offset_1 + a*een_rd_offset_2 + m*een_rd_offset_3 + nw*een_rd_offset_4]) * cn;
              }
            }

            cn = cn + cn;
            for(int j = 1; j<elec_num; j++) {
              factor_een_deriv_e[j + 3 + nw] = factor_een_deriv_e[j + 4 + nw] + (dtmp_c[j + 0*dtmp_offset_1 + a*dtmp_offset_2 + m*dtmp_offset_3 + k*dtmp_offset_4 + nw*dtmp_offset_5] * een_rescaled_n_deriv_e[j + 0 + a + (m+l) + nw] + dtmp_c[j + 1*dtmp_offset_1 + a*dtmp_offset_2 + m*dtmp_offset_3 + k*dtmp_offset_4 + nw*dtmp_offset_5] * een_rescaled_n_deriv_e[j + 1 + a + (m+l) + nw] + dtmp_c[j + 2*dtmp_offset_1 + a*dtmp_offset_2 + m*dtmp_offset_3 + k*dtmp_offset_4 + nw*dtmp_offset_5] * een_rescaled_n_deriv_e[j + 2 + a + (m+l) + nw] + dtmp_c[j + 0*dtmp_offset_1 + a*dtmp_offset_2 + (m+l)*dtmp_offset_3 + k*dtmp_offset_4 + nw*dtmp_offset_5] * een_rescaled_n_deriv_e[j + 0 + a + m + nw] + dtmp_c[j + 1*dtmp_offset_1 + a*dtmp_offset_2 + (m+l)*dtmp_offset_3 + k*dtmp_offset_4 + nw*dtmp_offset_5] * een_rescaled_n_deriv_e[j + 1 + a + m + nw] + dtmp_c[j + 2*dtmp_offset_1 + a*dtmp_offset_2 + (m+l)*dtmp_offset_3 + k*dtmp_offset_4 + nw*dtmp_offset_5] * een_rescaled_n_deriv_e[j + 2 + a + m + nw]) * cn;
            }
          }

        }
      }

      }

      return rc;
      }
#endif
