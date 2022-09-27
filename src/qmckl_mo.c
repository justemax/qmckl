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

#include "qmckl.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_ao_private_type.h"
#include "qmckl_ao_private_func.h"
#include "qmckl_mo_private_type.h"
#include "qmckl_mo_private_func.h"

qmckl_exit_code qmckl_init_mo_basis(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  ctx->mo_basis.uninitialized = (1 << 2) - 1;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_mo_basis_mo_num(qmckl_context context, const int64_t mo_num) {

  int32_t mask = 1 ;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_mo_*",
                             NULL);
   }

  if (mo_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_mo_basis_mo_num",
                           "mo_num <= 0");
  }

  ctx->mo_basis.mo_num = mo_num;

  ctx->mo_basis.uninitialized &= ~mask;
  ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
  if (ctx->mo_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_mo_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
  }
  return QMCKL_SUCCESS;
}

qmckl_exit_code  qmckl_set_mo_basis_coefficient(qmckl_context context, const double* coefficient) {

  int32_t mask = 1 << 1;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_mo_*",
                             NULL);
   }

  if (ctx->mo_basis.coefficient != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->mo_basis.coefficient);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_mo_basis_coefficient",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc(context, mem_info);
  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_mo_basis_coefficient",
                           NULL);
  }

  memcpy(new_array, coefficient, mem_info.size);

  ctx->mo_basis.coefficient = new_array;

  ctx->mo_basis.uninitialized &= ~mask;
  ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
  if (ctx->mo_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_mo_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
  }
  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code  qmckl_set_mo_basis_coefficient_device (qmckl_context context, const double  * coefficient, int device_id)
{

  int32_t mask = 1 << 1;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_mo_*",
                             NULL);
   }

  if (ctx->mo_basis.coefficient_device != NULL) {
    qmckl_exit_code rc = qmckl_free_device(context, ctx->mo_basis.coefficient_device, device_id);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_mo_basis_coefficient_device",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc_device(context, mem_info, device_id);
  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_mo_basis_coefficient_device",
                           NULL);
  }

  omp_target_memcpy(
    new_array, coefficient,
    mem_info.size,
    0, 0,
    device_id, omp_get_initial_device()
  );

  ctx->mo_basis.coefficient_device = new_array;

  ctx->mo_basis.uninitialized &= ~mask;
  ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
  if (ctx->mo_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_mo_basis_device(context, device_id);
    if (rc_ != QMCKL_SUCCESS) return rc_;
  }
  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code qmckl_finalize_mo_basis(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_finalize_mo_basis",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->mo_basis.coefficient_t != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->mo_basis.coefficient_t);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_finalize_mo_basis",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc(context, mem_info);
  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_finalize_mo_basis",
                           NULL);
  }

  assert (ctx->mo_basis.coefficient != NULL);

  for (int64_t i=0 ; i<ctx->ao_basis.ao_num ; ++i) {
    for (int64_t j=0 ; j<ctx->mo_basis.mo_num ; ++j) {
      new_array[i*ctx->mo_basis.mo_num + j] = ctx->mo_basis.coefficient[j*ctx->ao_basis.ao_num + i];
    }
  }

  ctx->mo_basis.coefficient_t = new_array;

  qmckl_exit_code rc; 
  if (ctx->mo_basis.mo_vgl != NULL) {
    rc = qmckl_free(context, ctx->mo_basis.mo_vgl);
    if (rc != QMCKL_SUCCESS) return rc;
    ctx->mo_basis.mo_vgl = NULL;
    ctx->mo_basis.mo_vgl_date = 0;
  }

  if (ctx->mo_basis.mo_value != NULL) {
    rc = qmckl_free(context, ctx->mo_basis.mo_value);
    if (rc != QMCKL_SUCCESS) return rc;
    ctx->mo_basis.mo_value = NULL;
    ctx->mo_basis.mo_value_date = 0;
  }

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_finalize_mo_basis_device(qmckl_context context, int device_id) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_finalize_mo_basis",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc_device(context, mem_info, device_id);
  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_finalize_mo_basis",
                           NULL);
  }

  assert (ctx->mo_basis.coefficient_device != NULL);

  if (ctx->mo_basis.coefficient_t_device != NULL) {
    qmckl_exit_code rc = qmckl_free_device(context, ctx->mo_basis.coefficient_device, device_id);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_finalize_mo_basis",
                             NULL);
    }
  }

  double * coefficient = ctx->mo_basis.coefficient_device;

  int64_t ao_num = ctx->ao_basis.ao_num;
  int64_t mo_num = ctx->mo_basis.mo_num;

  #pragma omp target is_device_ptr(new_array, coefficient)
  {
  #pragma omp parallel for collapse(2)
  for (int64_t i=0 ; i<ao_num ; ++i) {
    for (int64_t j=0 ; j<mo_num ; ++j) {
      new_array[i*mo_num + j] = coefficient[j*ao_num + i];
    }
  }
  }

  ctx->mo_basis.coefficient_t_device = new_array;
  qmckl_exit_code rc = QMCKL_SUCCESS;
  return rc;
}
#endif

qmckl_exit_code
qmckl_get_mo_basis_mo_num (const qmckl_context context,
                           int64_t* mo_num)
{
   if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_mo_num",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1;

  if ( (ctx->mo_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_mo_basis_mo_num",
                           NULL);
  }

  assert (ctx->mo_basis.mo_num > (int64_t) 0);
  *mo_num = ctx->mo_basis.mo_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_mo_basis_coefficient (const qmckl_context context,
                                double* const coefficient,
                                const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_coefficient",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 1;

  if ( (ctx->mo_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_mo_basis_coefficient",
                           NULL);
  }

  if (coefficient == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_mo_basis_coefficient",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.ao_num * ctx->mo_basis.mo_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_coefficient",
                           "Array too small: expected mo_num * ao_num.");
  }

  assert (ctx->mo_basis.coefficient != NULL);
  memcpy(coefficient, ctx->mo_basis.coefficient,
         ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double));

  return QMCKL_SUCCESS;
}

bool qmckl_mo_basis_provided(const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  return ctx->mo_basis.provided;
}

bool qmckl_mo_basis_select_mo (const qmckl_context context,
                               const int32_t* keep,
                               const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_select_mo",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if ( !(qmckl_mo_basis_provided(context)) ) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_mo_basis_select_mo",
                           NULL);
  }

  if (keep == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_mo_basis_select_mo",
                           "NULL pointer");
  }

  const int64_t mo_num = ctx->mo_basis.mo_num; 
  const int64_t ao_num = ctx->ao_basis.ao_num; 

  if (size_max < mo_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_select_mo",
                           "Array too small: expected mo_num.");
  }

  int64_t mo_num_new = 0;
  for (int64_t i=0 ; i<mo_num ; ++i) {
    if (keep[i] != 0) ++mo_num_new;
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ao_num * mo_num_new * sizeof(double);
  double* restrict coefficient   = (double*) qmckl_malloc(context, mem_info);

  int64_t k = 0;
  for (int64_t i=0 ; i<mo_num ; ++i) {
    if (keep[i] != 0) {
      memcpy(&(coefficient[k*ao_num]), &(ctx->mo_basis.coefficient[i*ao_num]), ao_num*sizeof(double));
      ++k;
    }
  }

  qmckl_exit_code rc = qmckl_free(context, ctx->mo_basis.coefficient);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->mo_basis.coefficient = coefficient;
  ctx->mo_basis.mo_num = mo_num_new;

  rc = qmckl_finalize_mo_basis(context);
  return rc;
}

qmckl_exit_code
qmckl_get_mo_basis_mo_value(qmckl_context context,
                            double* const mo_value,
                            const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_mo_basis_mo_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->point.num * ctx->mo_basis.mo_num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_mo_value",
                           "input array too small");
  }
  memcpy(mo_value, ctx->mo_basis.mo_value, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_mo_basis_mo_value_inplace (qmckl_context context,
                                     double* const mo_value,
                                     const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_mo_value",
                           NULL);
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->mo_basis.mo_num * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_mo_value",
                           "input array too small");
  }

  rc = qmckl_context_touch(context);
  if (rc != QMCKL_SUCCESS) return rc;

  double* old_array = ctx->mo_basis.mo_value;

  ctx->mo_basis.mo_value = mo_value;

  rc = qmckl_provide_mo_basis_mo_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->mo_basis.mo_value = old_array;

  return QMCKL_SUCCESS;
}



/* #+CALL: write_provider_pre( group="mo_basis", data="mo_value", dimension="ctx->mo_basis.mo_num * ctx->point.num") */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_value(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_mo_basis_mo_value",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_mo_basis_mo_value",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->mo_basis.mo_value_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->mo_basis.mo_num * ctx->point.num * sizeof(double);

    if (ctx->mo_basis.mo_value != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->mo_basis.mo_value, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->mo_basis.mo_value);
        assert (rc == QMCKL_SUCCESS);
        ctx->mo_basis.mo_value = NULL;
      }
    }

    /* Allocate array */
    if (ctx->mo_basis.mo_value == NULL) {

      double* mo_value = (double*) qmckl_malloc(context, mem_info);

      if (mo_value == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_mo_basis_mo_value",
                               NULL);
      }
      ctx->mo_basis.mo_value = mo_value;
    }

if (ctx->mo_basis.mo_vgl_date == ctx->point.date) {

  // mo_vgl has been computed at this step: Just copy the data.

  double * v = &(ctx->mo_basis.mo_value[0]);
  double * vgl = &(ctx->mo_basis.mo_vgl[0]);
  for (int i=0 ; i<ctx->point.num ; ++i) {
    for (int k=0 ; k<ctx->mo_basis.mo_num ; ++k) {
      v[k] = vgl[k];
    }
    v   += ctx->mo_basis.mo_num;
    vgl += ctx->mo_basis.mo_num * 5;
  }

} else {

  rc = qmckl_provide_ao_basis_ao_value(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_ao_value",
                           NULL);
  }

  rc = qmckl_compute_mo_basis_mo_value(context,
                                       ctx->ao_basis.ao_num,
                                       ctx->mo_basis.mo_num,
                                       ctx->point.num,
                                       ctx->mo_basis.coefficient_t,
                                       ctx->ao_basis.ao_value,
                                       ctx->mo_basis.mo_value);

}





/* #+CALL: write_provider_post( group="mo_basis", data="mo_value" ) */

/* #+RESULTS: */

if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->mo_basis.mo_value_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_mo_basis_mo_value (const qmckl_context context,
                                 const int64_t ao_num,
                                 const int64_t mo_num,
                                 const int64_t point_num,
                                 const double* coefficient_t,
                                 const double* ao_value,
                                 double* const mo_value )
{
#ifdef HAVE_HPC
  return qmckl_compute_mo_basis_mo_value_hpc (context, ao_num, mo_num, point_num, coefficient_t, ao_value, mo_value);
#else
  return qmckl_compute_mo_basis_mo_value_doc (context, ao_num, mo_num, point_num, coefficient_t, ao_value, mo_value);
#endif
}

#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_value_hpc (const qmckl_context context,
                                     const int64_t ao_num,
                                     const int64_t mo_num,
                                     const int64_t point_num,
                                     const double* restrict coefficient_t,
                                     const double* restrict ao_value,
                                     double* restrict const mo_value )
{
  assert (context != QMCKL_NULL_CONTEXT);

#ifdef HAVE_OPENMP
  #pragma omp parallel for
#endif
  for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {
    double* restrict const vgl1  = &(mo_value[ipoint*mo_num]);
    const double* restrict avgl1 = &(ao_value[ipoint*ao_num]);

    for (int64_t i=0 ; i<mo_num ; ++i) {
      vgl1[i] = 0.;
    }

    int64_t nidx=0;
    int64_t idx[ao_num];
    double  av1[ao_num];
    for (int64_t k=0 ; k<ao_num ; ++k) {
      if (avgl1[k] != 0.) {
        idx[nidx] = k;
        av1[nidx] = avgl1[k];
        ++nidx;
      }
    }

    int64_t n=0;

    for (n=0 ; n < nidx-4 ; n+=4) {
      const double* restrict ck1 = coefficient_t + idx[n  ]*mo_num;
      const double* restrict ck2 = coefficient_t + idx[n+1]*mo_num;
      const double* restrict ck3 = coefficient_t + idx[n+2]*mo_num;
      const double* restrict ck4 = coefficient_t + idx[n+3]*mo_num;

      const double a11 = av1[n  ];
      const double a21 = av1[n+1];
      const double a31 = av1[n+2];
      const double a41 = av1[n+3];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = vgl1[i] + ck1[i] * a11 + ck2[i] * a21 + ck3[i] * a31 + ck4[i] * a41;
      }
    }

    for (int64_t m=n ; m < nidx ; m+=1) {
      const double* restrict ck = coefficient_t + idx[m]*mo_num;
      const double a1 = av1[m];

#ifdef HAVE_OPENMP
  #pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] += ck[i] * a1;
      }
    }
  }
  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code
qmckl_get_mo_basis_mo_vgl(qmckl_context context,
                          double* const mo_vgl,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->point.num * 5 * ctx->mo_basis.mo_num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_mo_vgl",
                           "input array too small");
  }
  memcpy(mo_vgl, ctx->mo_basis.mo_vgl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_mo_basis_mo_vgl_inplace (qmckl_context context,
                                   double* const mo_vgl,
                                   const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_mo_vgl",
                           NULL);
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->mo_basis.mo_num * 5 * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_mo_vgl",
                           "input array too small");
  }

  rc = qmckl_context_touch(context);
  if (rc != QMCKL_SUCCESS) return rc;

  double* old_array = ctx->mo_basis.mo_vgl;

  ctx->mo_basis.mo_vgl = mo_vgl;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->mo_basis.mo_vgl = old_array;

  return QMCKL_SUCCESS;
}



/* #+CALL: write_provider_pre( group="mo_basis", data="mo_vgl", dimension="5 * ctx->mo_basis.mo_num * ctx->point.num") */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_vgl(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_mo_basis_mo_vgl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_mo_basis_mo_vgl",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->mo_basis.mo_vgl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 5 * ctx->mo_basis.mo_num * ctx->point.num * sizeof(double);

    if (ctx->mo_basis.mo_vgl != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->mo_basis.mo_vgl, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->mo_basis.mo_vgl);
        assert (rc == QMCKL_SUCCESS);
        ctx->mo_basis.mo_vgl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->mo_basis.mo_vgl == NULL) {

      double* mo_vgl = (double*) qmckl_malloc(context, mem_info);

      if (mo_vgl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_mo_basis_mo_vgl",
                               NULL);
      }
      ctx->mo_basis.mo_vgl = mo_vgl;
    }

rc = qmckl_provide_ao_basis_ao_vgl(context);
if (rc != QMCKL_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_NOT_PROVIDED,
                         "qmckl_ao_basis",
                         NULL);
}

rc = qmckl_compute_mo_basis_mo_vgl(context,
                                   ctx->ao_basis.ao_num,
                                   ctx->mo_basis.mo_num,
                                   ctx->point.num,
                                   ctx->mo_basis.coefficient_t,
                                   ctx->ao_basis.ao_vgl,
                                   ctx->mo_basis.mo_vgl);



/* #+CALL: write_provider_post( group="mo_basis", data="mo_vgl" ) */

/* #+RESULTS: */

if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->mo_basis.mo_vgl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl (const qmckl_context context,
                            const int64_t ao_num,
                            const int64_t mo_num,
                            const int64_t point_num,
                            const double* coefficient_t,
                            const double* ao_vgl,
                            double* const mo_vgl )
{
#ifdef HAVE_HPC
  return qmckl_compute_mo_basis_mo_vgl_hpc (context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl);
#else
  return qmckl_compute_mo_basis_mo_vgl_doc (context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl);
#endif
}

#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl_hpc (const qmckl_context context,
                                   const int64_t ao_num,
                                   const int64_t mo_num,
                                   const int64_t point_num,
                                   const double* restrict coefficient_t,
                                   const double* restrict ao_vgl,
                                   double* restrict const mo_vgl )
{
  assert (context != QMCKL_NULL_CONTEXT);

#ifdef HAVE_OPENMP
  #pragma omp parallel for
#endif
  for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {
    double* restrict const vgl1 = &(mo_vgl[ipoint*5*mo_num]);
    double* restrict const vgl2 =  vgl1 + mo_num;
    double* restrict const vgl3 =  vgl1 + (mo_num << 1);
    double* restrict const vgl4 =  vgl1 + (mo_num << 1) + mo_num;
    double* restrict const vgl5 =  vgl1 + (mo_num << 2);

    const double* restrict avgl1 = &(ao_vgl[ipoint*5*ao_num]);
    const double* restrict avgl2 = avgl1 + ao_num;
    const double* restrict avgl3 = avgl1 + (ao_num << 1);
    const double* restrict avgl4 = avgl1 + (ao_num << 1) + ao_num;
    const double* restrict avgl5 = avgl1 + (ao_num << 2);

    for (int64_t i=0 ; i<mo_num ; ++i) {
      vgl1[i] = 0.;
      vgl2[i] = 0.;
      vgl3[i] = 0.;
      vgl4[i] = 0.;
      vgl5[i] = 0.;
    }

    int64_t nidx=0;
    int64_t idx[ao_num];
    double  av1[ao_num];
    double  av2[ao_num];
    double  av3[ao_num];
    double  av4[ao_num];
    double  av5[ao_num];
    for (int64_t k=0 ; k<ao_num ; ++k) {
      if (avgl1[k] != 0.) {
        idx[nidx] = k;
        av1[nidx] = avgl1[k];
        av2[nidx] = avgl2[k];
        av3[nidx] = avgl3[k];
        av4[nidx] = avgl4[k];
        av5[nidx] = avgl5[k];
        ++nidx;
      }
    }

    int64_t n=0;

    for (n=0 ; n < nidx-4 ; n+=4) {
      const double* restrict ck1 = coefficient_t + idx[n  ]*mo_num;
      const double* restrict ck2 = coefficient_t + idx[n+1]*mo_num;
      const double* restrict ck3 = coefficient_t + idx[n+2]*mo_num;
      const double* restrict ck4 = coefficient_t + idx[n+3]*mo_num;

      const double a11 = av1[n  ];
      const double a21 = av1[n+1];
      const double a31 = av1[n+2];
      const double a41 = av1[n+3];

      const double a12 = av2[n  ];
      const double a22 = av2[n+1];
      const double a32 = av2[n+2];
      const double a42 = av2[n+3];

      const double a13 = av3[n  ];
      const double a23 = av3[n+1];
      const double a33 = av3[n+2];
      const double a43 = av3[n+3];

      const double a14 = av4[n  ];
      const double a24 = av4[n+1];
      const double a34 = av4[n+2];
      const double a44 = av4[n+3];

      const double a15 = av5[n  ];
      const double a25 = av5[n+1];
      const double a35 = av5[n+2];
      const double a45 = av5[n+3];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = vgl1[i] + ck1[i] * a11 + ck2[i] * a21 + ck3[i] * a31 + ck4[i] * a41;
        vgl2[i] = vgl2[i] + ck1[i] * a12 + ck2[i] * a22 + ck3[i] * a32 + ck4[i] * a42;
        vgl3[i] = vgl3[i] + ck1[i] * a13 + ck2[i] * a23 + ck3[i] * a33 + ck4[i] * a43;
        vgl4[i] = vgl4[i] + ck1[i] * a14 + ck2[i] * a24 + ck3[i] * a34 + ck4[i] * a44;
        vgl5[i] = vgl5[i] + ck1[i] * a15 + ck2[i] * a25 + ck3[i] * a35 + ck4[i] * a45;
      }
    }

    for (int64_t m=n ; m < nidx ; m+=1) {
      const double* restrict ck = coefficient_t + idx[m]*mo_num;
      const double a1 = av1[m];
      const double a2 = av2[m];
      const double a3 = av3[m];
      const double a4 = av4[m];
      const double a5 = av5[m];

#ifdef HAVE_OPENMP
  #pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] += ck[i] * a1;
        vgl2[i] += ck[i] * a2;
        vgl3[i] += ck[i] * a3;
        vgl4[i] += ck[i] * a4;
        vgl5[i] += ck[i] * a5;
      }
    }
  }
  return QMCKL_SUCCESS;
}
#endif
