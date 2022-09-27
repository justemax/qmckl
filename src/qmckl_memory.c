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
#include <assert.h>

#include <omp.h>

#include "qmckl.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_func.h"

void* qmckl_malloc(qmckl_context context, const qmckl_memory_info_struct info) {

  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  /* Allocate memory and zero it */
#ifdef HAVE_HPC
  void * pointer = aligned_alloc(64, ((info.size+64) >> 6) << 6 );
#else
  void * pointer = malloc(info.size);
#endif
  if (pointer == NULL) {
    return NULL;
  }
  memset(pointer, 0, info.size);

  qmckl_lock(context);
  {
    /* If qmckl_memory_struct is full, reallocate a larger one */
    if (ctx->memory.n_allocated == ctx->memory.array_size) {
      const size_t old_size = ctx->memory.array_size;
      qmckl_memory_info_struct * new_array = realloc(ctx->memory.element,
                                                          2L * old_size *
                                                          sizeof(qmckl_memory_info_struct));
      if (new_array == NULL) {
        qmckl_unlock(context);
        free(pointer);
        return NULL;
      }

      memset( &(new_array[old_size]), 0, old_size * sizeof(qmckl_memory_info_struct) );
      ctx->memory.element = new_array;
      ctx->memory.array_size = 2L * old_size;
    }

    /* Find first NULL entry */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].size > (size_t) 0) {
      pos += (size_t) 1;
    }
    assert (ctx->memory.element[pos].size == (size_t) 0);

    /* Copy info at the new location */
    memcpy(&(ctx->memory.element[pos]), &info, sizeof(qmckl_memory_info_struct));
    ctx->memory.element[pos].pointer = pointer;
    ctx->memory.n_allocated += (size_t) 1;
  }
  qmckl_unlock(context);

  return pointer;
}

#ifdef HAVE_DEVICE_POINTERS
void* qmckl_malloc_device(qmckl_context context, const qmckl_memory_info_struct info, int device_id) {

  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  /* Allocate memory and zero it */
  void * pointer = omp_target_alloc(info.size, device_id);
  if (pointer == NULL) {
    return NULL;
  }

  // TODO
  // Memset to 0 of size info.size
  // memset(pointer, 0, info.size);


  qmckl_lock(context);
  {
    /* If qmckl_memory_struct is full, reallocate a larger one */
    if (ctx->memory.n_allocated == ctx->memory.array_size) {
      const size_t old_size = ctx->memory.array_size;
      qmckl_memory_info_struct * new_array = realloc(ctx->memory.element,
                                                          2L * old_size *
                                                          sizeof(qmckl_memory_info_struct));
      if (new_array == NULL) {
        qmckl_unlock(context);
        free(pointer);
        return NULL;
      }

      memset( &(new_array[old_size]), 0, old_size * sizeof(qmckl_memory_info_struct) );
      ctx->memory.element = new_array;
      ctx->memory.array_size = 2L * old_size;
    }

    /* Find first NULL entry */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].size > (size_t) 0) {
      pos += (size_t) 1;
    }
    assert (ctx->memory.element[pos].size == (size_t) 0);

    /* Copy info at the new location */
    ctx->memory.element[pos].size    = info.size;
    ctx->memory.element[pos].pointer = pointer;
    ctx->memory.n_allocated += (size_t) 1;
  }
  qmckl_unlock(context);

  return pointer;
}
#endif

qmckl_exit_code qmckl_free(qmckl_context context, void * const ptr) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith(context,
                            QMCKL_INVALID_CONTEXT,
                            "qmckl_free",
                            NULL);
  }

  if (ptr == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_free",
                          "NULL pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  qmckl_lock(context);
  {
    /* Find pointer in array of saved pointers */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].pointer != ptr) {
      pos += (size_t) 1;
    }

    if (pos >= ctx->memory.array_size) {
      /* Not found */
      qmckl_unlock(context);
      return qmckl_failwith(context,
                            QMCKL_INVALID_ARG_2,
                            "qmckl_free",
                            "Pointer not found in context");
    }

    free(ptr);

    memset( &(ctx->memory.element[pos]), 0, sizeof(qmckl_memory_info_struct) );
    ctx->memory.n_allocated -= (size_t) 1;
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code qmckl_free_device(qmckl_context context, void * const ptr, int device_id) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith(context,
                            QMCKL_INVALID_CONTEXT,
                            "qmckl_free",
                            NULL);
  }

  if (ptr == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_free",
                          "NULL pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  qmckl_lock(context);
  {
    /* Find pointer in array of saved pointers */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].pointer != ptr) {
      pos += (size_t) 1;
    }

    if (pos >= ctx->memory.array_size) {
      /* Not found */
      qmckl_unlock(context);
      return qmckl_failwith(context,
                            QMCKL_FAILURE,
                            "qmckl_free",
                            "Pointer not found in context");
    }

    omp_target_free(ptr, device_id);

    memset( &(ctx->memory.element[pos]), 0, sizeof(qmckl_memory_info_struct) );
    ctx->memory.n_allocated -= (size_t) 1;
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code
qmckl_get_malloc_info(qmckl_context context,
                      const void* ptr, 
                      qmckl_memory_info_struct* info) 
{

  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  if (ptr == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_get_malloc_info",
                          "Null pointer");
  }

  if (info == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_3,
                          "qmckl_get_malloc_info",
                          "Null pointer");
  }

  qmckl_lock(context);
  {
    /* Find the pointer entry */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].pointer != ptr) {
      pos += (size_t) 1;
    }

    if (pos >= ctx->memory.array_size) {
      /* Not found */
      qmckl_unlock(context);
      return qmckl_failwith(context,
                            QMCKL_INVALID_ARG_2,
                            "qmckl_get_malloc_info",
                            "Pointer not found in context");
    }

    /* Copy info */
    memcpy(info, &(ctx->memory.element[pos]), sizeof(qmckl_memory_info_struct));
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}
