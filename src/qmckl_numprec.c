#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "qmckl.h"
#include "qmckl_context_private_type.h"

qmckl_exit_code qmckl_set_numprec_precision(const qmckl_context context, const int precision) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT)
    return QMCKL_INVALID_CONTEXT;

  if (precision <  2) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_update_numprec_precision",
                          "precision < 2");
  }

  if (precision > 53) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_update_numprec_precision",
                          "precision > 53");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  /* This should be always true because the context is valid */
  assert (ctx != NULL);

  qmckl_lock(context);
  {
    ctx->numprec.precision = (uint32_t) precision;
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}

int qmckl_get_numprec_precision(const qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith(context,
                      QMCKL_INVALID_CONTEXT,
                      "qmckl_get_numprec_precision",
                      "");
  }

  const qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  return ctx->numprec.precision;
}

qmckl_exit_code qmckl_set_numprec_range(const qmckl_context context, const int range) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT)
    return QMCKL_INVALID_CONTEXT;

  if (range <  2) {
    return qmckl_failwith(context,
                    QMCKL_INVALID_ARG_2,
                    "qmckl_set_numprec_range",
                    "range < 2");
  }

  if (range > 11) {
    return qmckl_failwith(context,
                    QMCKL_INVALID_ARG_2,
                    "qmckl_set_numprec_range",
                    "range > 11");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  /* This should be always true because the context is valid */
  assert (ctx != NULL);

  qmckl_lock(context);
  {
    ctx->numprec.range = (uint32_t) range;
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}

int qmckl_get_numprec_range(const qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith(context,
                      QMCKL_INVALID_CONTEXT,
                      "qmckl_get_numprec_range",
                      "");
  }

  const qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  return ctx->numprec.range;
}

double qmckl_get_numprec_epsilon(const qmckl_context context) {
  const int precision = qmckl_get_numprec_precision(context);
  return 1. /  (double) (1L << (precision-2));
}
