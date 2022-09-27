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
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>

#include "qmckl.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_func.h"

qmckl_context qmckl_context_check(const qmckl_context context) {

  if (context == QMCKL_NULL_CONTEXT)
    return QMCKL_NULL_CONTEXT;

  const qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  /* Try to access memory */
  if (ctx->tag != VALID_TAG) {
      return QMCKL_NULL_CONTEXT;
  }

  return context;
}

qmckl_exit_code
qmckl_context_touch(const qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_context_touch",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  ctx->date += 1UL;
  ctx->point.date = ctx-> date;
  return QMCKL_SUCCESS;
}

qmckl_context qmckl_context_create() {

  qmckl_context_struct* const ctx =
    (qmckl_context_struct*) malloc (sizeof(qmckl_context_struct));

  if (ctx == NULL) {
    return QMCKL_NULL_CONTEXT;
  }

  /* Set all pointers and values to NULL */
  {
    memset(ctx, 0, sizeof(qmckl_context_struct));
  }

  /* Initialize lock */
  {
    pthread_mutexattr_t attr;
    int rc;

    rc = pthread_mutexattr_init(&attr);
    assert (rc == 0);

#ifdef PTHREAD_MUTEX_RECURSIVE
    (void) pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
#endif

    rc = pthread_mutex_init ( &(ctx->mutex), &attr);
    assert (rc == 0);

    (void) pthread_mutexattr_destroy(&attr);
  }

  /* Initialize data */
  {
    ctx->tag = VALID_TAG;

    const qmckl_context context = (qmckl_context) ctx;
    assert ( qmckl_context_check(context) != QMCKL_NULL_CONTEXT );

    qmckl_exit_code rc;

    ctx->numprec.precision = QMCKL_DEFAULT_PRECISION;
    ctx->numprec.range = QMCKL_DEFAULT_RANGE;

    rc = qmckl_init_point(context);
    assert (rc == QMCKL_SUCCESS);

    rc = qmckl_init_electron(context);
    assert (rc == QMCKL_SUCCESS);

    rc = qmckl_init_nucleus(context);
    assert (rc == QMCKL_SUCCESS);

    rc = qmckl_init_ao_basis(context);
    assert (rc == QMCKL_SUCCESS);

    rc = qmckl_init_mo_basis(context);
    assert (rc == QMCKL_SUCCESS);

    rc = qmckl_init_determinant(context);
    assert (rc == QMCKL_SUCCESS);
  }

  /* Allocate qmckl_memory_struct */
  {
    const size_t size = 128L;
    qmckl_memory_info_struct * new_array = calloc(size, sizeof(qmckl_memory_info_struct));
    if (new_array == NULL) {
      free(ctx);
      return QMCKL_NULL_CONTEXT;
    }
    memset( &(new_array[0]), 0, size * sizeof(qmckl_memory_info_struct) );

    ctx->memory.element = new_array;
    ctx->memory.array_size = size;
    ctx->memory.n_allocated = (size_t) 0;
  }

  return (qmckl_context) ctx;
}

void qmckl_lock(qmckl_context context) {
  if (context == QMCKL_NULL_CONTEXT)
    return ;
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  errno = 0;
  int rc = pthread_mutex_lock( &(ctx->mutex) );
  if (rc != 0) {
    fprintf(stderr, "DEBUG qmckl_lock:%s\n", strerror(rc) );
    fflush(stderr);
  }
  assert (rc == 0);
  ctx->lock_count += 1;
}

void qmckl_unlock(const qmckl_context context) {
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  int rc = pthread_mutex_unlock( &(ctx->mutex) );
  if (rc != 0) {
    fprintf(stderr, "DEBUG qmckl_unlock:%s\n", strerror(rc) );
    fflush(stderr);
  }
  assert (rc == 0);
  ctx->lock_count -= 1;
}

qmckl_context qmckl_context_copy(const qmckl_context context) {

  const qmckl_context checked_context = qmckl_context_check(context);

  if (checked_context == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  /*
  qmckl_lock(context);
  {

    const qmckl_context_struct* const old_ctx =
      (qmckl_context_struct*) checked_context;

    qmckl_context_struct* const new_ctx =
      (qmckl_context_struct*) malloc (context, sizeof(qmckl_context_struct));

    if (new_ctx == NULL) {
      qmckl_unlock(context);
      return QMCKL_NULL_CONTEXT;
    }

     * Copy the old context on the new one *
     * TODO Deep copies should be done here *
    memcpy(new_ctx, old_ctx, sizeof(qmckl_context_struct));

    qmckl_unlock( (qmckl_context) new_ctx );

    return (qmckl_context) new_ctx;
  }
  qmckl_unlock(context);
*/
    return QMCKL_NULL_CONTEXT;
}

qmckl_exit_code
qmckl_context_destroy (const qmckl_context context)
{

  const qmckl_context checked_context = qmckl_context_check(context);
  if (checked_context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);  /* Shouldn't be possible because the context is valid */

  qmckl_lock(context);
  {
    /* Memory: Remove all allocated data */
    for (size_t pos = (size_t) 0 ; pos < ctx->memory.array_size ; ++pos) {
      if (ctx->memory.element[pos].pointer != NULL) {
        free(ctx->memory.element[pos].pointer);
        memset( &(ctx->memory.element[pos]), 0, sizeof(qmckl_memory_info_struct) );
        ctx->memory.n_allocated -= 1;
      }
    }
    assert (ctx->memory.n_allocated == (size_t) 0);
    free(ctx->memory.element);
    ctx->memory.element = NULL;
    ctx->memory.array_size = (size_t) 0;
  }
  qmckl_unlock(context);

  ctx->tag = INVALID_TAG;

  const int rc_destroy = pthread_mutex_destroy( &(ctx->mutex) );
  if (rc_destroy != 0) {
/* DEBUG */
     fprintf(stderr, "qmckl_context_destroy: %s (count = %d)\n", strerror(rc_destroy), ctx->lock_count);
     abort();
  }

  free(ctx);

  return QMCKL_SUCCESS;
}
