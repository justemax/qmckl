#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#ifdef HAVE_DEVICE_POINTERS
#include <omp.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "qmckl.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_blas_private_type.h"

#include "qmckl_memory_private_func.h"
#include "qmckl_blas_private_func.h"




/* Allocates a new vector. If the allocation failed the size is zero. */


qmckl_vector
qmckl_vector_alloc( qmckl_context context,
                    const int64_t size)
{
  /* Should always be true by contruction */
  assert (size > (int64_t) 0);

  qmckl_vector result;
  result.size = size;

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = size * sizeof(double);
  result.data = (double*) qmckl_malloc (context, mem_info);

  if (result.data == NULL) {
    result.size = (int64_t) 0;
  }

  return result;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_vector
qmckl_vector_alloc_device( qmckl_context context,
                           const int64_t size,
                           int device_id)
{
  /* Should always be true by contruction */
  assert (size > (int64_t) 0);

  qmckl_vector result;
  result.size = size;

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = size * sizeof(double);
  result.data_device = (double*) qmckl_malloc_device (context, mem_info, device_id);

  if (result.data == NULL) {
    result.size = (int64_t) 0;
  }

  return result;
}
#endif

qmckl_exit_code
qmckl_vector_free( qmckl_context context,
                   qmckl_vector* vector)
{
  if (vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_vector_free",
                           "Null pointer");
  }

  /* Always true */
  assert (vector->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free(context, vector->data);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  vector->size = (int64_t) 0;
  vector->data = NULL;
  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_vector_free_device( qmckl_context context,
                          qmckl_vector* vector,
                          int device_id)
{
  /* Always true */
  assert (vector->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free_device(context, vector->data_device, device_id);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  vector->size = (int64_t) 0;
  vector->data = NULL;
  return QMCKL_SUCCESS;
}
#endif




/* Allocates a new matrix. If the allocation failed the sizes are zero. */


qmckl_matrix
qmckl_matrix_alloc( qmckl_context context,
                    const int64_t size1,
                    const int64_t size2)
{
  /* Should always be true by contruction */
  assert (size1 * size2 > (int64_t) 0);

  qmckl_matrix result;

  result.size[0] = size1;
  result.size[1] = size2;

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = size1 * size2 * sizeof(double);
  result.data = (double*) qmckl_malloc (context, mem_info);

  if (result.data == NULL) {
    result.size[0] = (int64_t) 0;
    result.size[1] = (int64_t) 0;
  }

  return result;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_matrix
qmckl_matrix_alloc_device( qmckl_context context,
                           const int64_t size1,
                           const int64_t size2,
                           int device_id)
{
  /* Should always be true by contruction */
  assert (size1 * size2 > (int64_t) 0);

  qmckl_matrix result;

  result.size[0] = size1;
  result.size[1] = size2;

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = size1 * size2 * sizeof(double);
  result.data_device = (double*) qmckl_malloc_device (context, mem_info, device_id);

  if (result.data_device == NULL) {
    result.size[0] = (int64_t) 0;
    result.size[1] = (int64_t) 0;
  }

  return result;
}
#endif

qmckl_exit_code
qmckl_matrix_free( qmckl_context context,
                   qmckl_matrix* matrix)
{
  if (matrix == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_matrix_free",
                           "Null pointer");
  }

  /* Always true */
  assert (matrix->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free(context, matrix->data);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }
  matrix->data = NULL;
  matrix->size[0] = (int64_t) 0;
  matrix->size[1] = (int64_t) 0;

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_matrix_free_device( qmckl_context context,
                          qmckl_matrix* matrix,
                          int device_id)
{
  /* Always true */
  assert (matrix->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free_device(context, matrix->data, device_id);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }
  matrix->data = NULL;
  matrix->size[0] = (int64_t) 0;
  matrix->size[1] = (int64_t) 0;

  return QMCKL_SUCCESS;
}
#endif



/* Allocates memory for a tensor. If the allocation failed, the size */
/* is zero. */


qmckl_tensor
qmckl_tensor_alloc( qmckl_context context,
                    const int64_t  order,
                    const int64_t* size)
{
  /* Should always be true by contruction */
  assert (order > 0);
  assert (order <= QMCKL_TENSOR_ORDER_MAX);
  assert (size  != NULL);

  qmckl_tensor result;
  result.order = order;

  int64_t prod_size = (int64_t) 1;
  for (int64_t i=0 ; i<order ; ++i) {
    assert (size[i] > (int64_t) 0);
    result.size[i] = size[i];
    prod_size *= size[i];
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = prod_size * sizeof(double);

  result.data = (double*) qmckl_malloc (context, mem_info);

  if (result.data == NULL) {
    memset(&result, 0, sizeof(qmckl_tensor));
  }

  return result;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_tensor
qmckl_tensor_alloc_device( qmckl_context context,
                           const int64_t  order,
                           const int64_t* size,
	                       int device_id)
{
  /* Should always be true by contruction */
  assert (order > 0);
  assert (order <= QMCKL_TENSOR_ORDER_MAX);
  assert (size  != NULL);

  qmckl_tensor result;
  result.order = order;

  int64_t prod_size = (int64_t) 1;
  for (int64_t i=0 ; i<order ; ++i) {
    assert (size[i] > (int64_t) 0);
    result.size[i] = size[i];
    prod_size *= size[i];
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = prod_size * sizeof(double);

  result.data_device = (double*) qmckl_malloc_device (context, mem_info, device_id);

  if (result.data_device == NULL) {
    memset(&result, 0, sizeof(qmckl_tensor));
  }

  return result;
}
#endif

qmckl_exit_code
qmckl_tensor_free( qmckl_context context,
                   qmckl_tensor* tensor)
{
  if (tensor == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_tensor_free",
                           "Null pointer");
  }

  /* Always true */
  assert (tensor->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free(context, tensor->data);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  memset(tensor, 0, sizeof(qmckl_tensor));

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_tensor_free_device( qmckl_context context,
                          qmckl_tensor* tensor,
	                      int device_id)
{
  /* Always true */
  assert (tensor->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free_device(context, tensor->data_device, device_id);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  memset(tensor, 0, sizeof(qmckl_tensor));

  return QMCKL_SUCCESS;
}
#endif



/* Reshapes a vector into a matrix. */


qmckl_matrix
qmckl_matrix_of_vector(const qmckl_vector vector,
                       const int64_t size1,
                       const int64_t size2)
{
  /* Always true */
  assert (size1 * size2 == vector.size);

  qmckl_matrix result;

  result.size[0] = size1;
  result.size[1] = size2;
  result.data    = vector.data;

  return result;
}



/* Reshapes a vector into a tensor. */


qmckl_tensor
qmckl_tensor_of_vector(const qmckl_vector vector,
                       const int64_t order,
                       const int64_t* size)
{
  qmckl_tensor result;

  int64_t prod_size = 1;
  for (int64_t i=0 ; i<order ; ++i) {
    result.size[i] = size[i];
    prod_size *= size[i];
  }
  assert (prod_size == vector.size);

  result.data = vector.data;

  return result;
}



/* Reshapes a matrix into a vector. */


qmckl_vector
qmckl_vector_of_matrix(const qmckl_matrix matrix)
{
  qmckl_vector result;

  result.size = matrix.size[0] * matrix.size[1];
  result.data = matrix.data;

  return result;
}



/* Reshapes a matrix into a tensor. */


qmckl_tensor
qmckl_tensor_of_matrix(const qmckl_matrix matrix,
                       const int64_t order,
                       const int64_t* size)
{
  qmckl_tensor result;

  int64_t prod_size = 1;
  for (int64_t i=0 ; i<order ; ++i) {
    result.size[i] = size[i];
    prod_size *= size[i];
  }
  assert (prod_size == matrix.size[0] * matrix.size[1]);

  result.data = matrix.data;

  return result;
}



/* Reshapes a tensor into a vector. */


qmckl_vector
qmckl_vector_of_tensor(const qmckl_tensor tensor)
{
  int64_t prod_size = (int64_t) tensor.size[0];
  for (int64_t i=1 ; i<tensor.order ; i++) {
    prod_size *= tensor.size[i];
  }

  qmckl_vector result;

  result.size = prod_size;
  result.data = tensor.data;

  return result;
}



/* Reshapes a tensor into a vector. */


qmckl_matrix
qmckl_matrix_of_tensor(const qmckl_tensor tensor,
                       const int64_t size1,
                       const int64_t size2)
{
  /* Always true */
  int64_t prod_size = (int64_t) 1;
  for (int64_t i=0 ; i<tensor.order ; i++) {
    prod_size *= tensor.size[i];
  }
  assert (prod_size == size1 * size2);

  qmckl_matrix result;

  result.size[0] = size1;
  result.size[1] = size2;
  result.data = tensor.data;

  return result;
}

qmckl_vector
qmckl_vector_set(qmckl_vector vector, double value)
{
  for (int64_t i=0 ; i<vector.size ; ++i) {
    qmckl_vec(vector, i) = value;
  }
  return vector;
}

qmckl_matrix
qmckl_matrix_set(qmckl_matrix matrix, double value)
{
  qmckl_vector vector = qmckl_vector_of_matrix(matrix);
  for (int64_t i=0 ; i<vector.size ; ++i) {
    qmckl_vec(vector, i) = value;
  }
  return qmckl_matrix_of_vector(vector, matrix.size[0], matrix.size[1]);
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_matrix
qmckl_matrix_set_device(qmckl_matrix matrix, double value)
{
	// Recompute array size
	int prod_size = matrix.size[0] * matrix.size[1];

	double * data_device = matrix.data_device;
	#pragma omp target is_device_ptr(data_device)
	{
	for(int i=0; i<prod_size; i++) {
    data_device[i] = value;
	}
	}
  return matrix;
}
#endif

qmckl_tensor
qmckl_tensor_set(qmckl_tensor tensor, double value)
{
  qmckl_vector vector = qmckl_vector_of_tensor(tensor);
  for (int64_t i=0 ; i<vector.size ; ++i) {
    qmckl_vec(vector, i) = value;
  }
  return qmckl_tensor_of_vector(vector, tensor.order, tensor.size);
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_tensor
qmckl_tensor_set_device(qmckl_tensor tensor, double value)
{
	// Recompute array size
	int prod_size = 1;

	for(int i=0; i<tensor.order; i++) {
	  prod_size *= tensor.size[i];
	}

	double * data_device = tensor.data_device;
	#pragma omp target is_device_ptr(data_device)
	{
	for(int i=0; i<prod_size; i++) {
    data_device[i] = value;
	}
	}
  return tensor;
}
#endif



/* Converts a vector to a ~double*~. */


qmckl_exit_code
qmckl_double_of_vector(const qmckl_context context,
                       const qmckl_vector vector,
                       double* const target,
                       const int64_t size_max)
{
  /* Always true by construction */
  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);
  assert (vector.size > (int64_t) 0);
  assert (target != NULL);
  assert (size_max > (int64_t) 0);
  assert (size_max >= vector.size);
  for (int64_t i=0 ; i<vector.size ; ++i) {
    target[i] = vector.data[i];
  }
  return QMCKL_SUCCESS;

}



/* Converts a matrix to a ~double*~. */


qmckl_exit_code
qmckl_double_of_matrix(const qmckl_context context,
                       const qmckl_matrix matrix,
                       double* const target,
                       const int64_t size_max)
{
  qmckl_vector vector = qmckl_vector_of_matrix(matrix);
  return qmckl_double_of_vector(context, vector, target, size_max);
}



/* Converts a tensor to a ~double*~. */


qmckl_exit_code
qmckl_double_of_tensor(const qmckl_context context,
                       const qmckl_tensor tensor,
                       double* const target,
                       const int64_t size_max)
{
  qmckl_vector vector = qmckl_vector_of_tensor(tensor);
  return qmckl_double_of_vector(context, vector, target, size_max);
}



/* Converts a ~double*~ to a vector. */


qmckl_exit_code
qmckl_vector_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_vector* vector_out)
{
  qmckl_vector vector = *vector_out;
  /* Always true by construction */
  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  if (vector.size == 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_4,
                           "qmckl_vector_of_double",
                           "Vector not allocated");
  }

  if (vector.size != size_max) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_4,
                             "qmckl_vector_of_double",
                             "Wrong vector size");
  }

  for (int64_t i=0 ; i<vector.size ; ++i) {
    vector.data[i] = target[i];
  }

  *vector_out = vector;
  return QMCKL_SUCCESS;

}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_vector_of_double_device(const qmckl_context context,
                              const double* target,
                              const int64_t size_max,
                              qmckl_vector* vector_out,
                              int device_id)
{

  // Accepts an host array an copies it in the device section of vector_out
  // (assuming the vector is already allocated)

  qmckl_vector vector = *vector_out;
  /* Always true by construction */
  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  if (vector.size == 0) {
    // This error is thrown
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_4,
                           "qmckl_vector_of_double",
                           "Vector not allocated");
  }

  if (vector.size != size_max) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_4,
                             "qmckl_vector_of_double",
                             "Wrong vector size");
  }

  omp_target_memcpy(
    vector.data_device, target,
    vector.size * sizeof(double),
    0, 0,
    device_id, omp_get_initial_device()
  );

  *vector_out = vector;
  return QMCKL_SUCCESS;

}
#endif



/* Converts a ~double*~ to a matrix. */


qmckl_exit_code
qmckl_matrix_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_matrix* matrix)
{
  qmckl_vector vector = qmckl_vector_of_matrix(*matrix);
  qmckl_exit_code rc =
    qmckl_vector_of_double(context, target, size_max, &vector);
  *matrix = qmckl_matrix_of_vector(vector, matrix->size[0], matrix->size[1]);
  return rc;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_matrix_of_double_device(const qmckl_context context,
                              const double* target,
                              const int64_t size_max,
                              qmckl_matrix* matrix_out,
                              int device_id)
{

  // Accepts an host array an copies it in the device section of matrix
  // (assuming the matrix is already allocated)

  qmckl_matrix matrix = *matrix_out;
  /* Always true by construction */
  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  if (matrix.size[0] * matrix.size[1] == 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_4,
                           "qmckl_matrix_of_double",
                           "Matrix not allocated");
  }

  if (matrix.size[0] * matrix.size[1] != size_max) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_4,
                             "qmckl_matrix_of_double",
                             "Wrong vector size");
  }

  omp_target_memcpy(
    matrix.data_device, target,
    size_max * sizeof(double),
    0, 0,
    device_id, omp_get_initial_device()
  );

  *matrix_out = matrix;
  return QMCKL_SUCCESS;
}
#endif



/* Converts a ~double*~ to a tensor. */


qmckl_exit_code
qmckl_tensor_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_tensor* tensor)
{
  qmckl_vector vector = qmckl_vector_of_tensor(*tensor);
  qmckl_exit_code rc =
    qmckl_vector_of_double(context, target, size_max, &vector);
  *tensor = qmckl_tensor_of_vector(vector, tensor->order, tensor->size);
  return rc;
}

qmckl_exit_code
qmckl_matmul (const qmckl_context context,
              const char TransA,
              const char TransB,
              const double alpha,
              const qmckl_matrix A,
              const qmckl_matrix B,
              const double beta,
              qmckl_matrix* const C )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (TransA != 'N' && TransA != 'T') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_matmul",
                           "TransA should be 'N' or 'T'");
  }

  if (TransB != 'N' && TransB != 'T') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_matmul",
                           "TransB should be 'N' or 'T'");
  }

  if (A.size[0] < 1) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_5,
                           "qmckl_matmul",
                           "Invalid size for A");
  }

  if (B.size[0] < 1) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_6,
                           "qmckl_matmul",
                           "Invalid size for B");
  }

  if (C == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_8,
                           "qmckl_matmul",
                           "Null pointer");
  }

  int t = 0;
  if (TransA == 'T') t +=1;
  if (TransB == 'T') t +=2;
  /*
    | t | TransA | TransB |
    +---+--------+--------+
    | 0 | N      | N      |
    | 1 | T      | N      |
    | 2 | N      | T      |
    | 3 | T      | T      |
  */

  switch (t) {
  case 0:
    if (A.size[1] != B.size[0]) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_matmul",
                             "A and B have incompatible dimensions");
    }
    C->size[0] = A.size[0];
    C->size[1] = B.size[1];
    rc = qmckl_dgemm (context, 'N', 'N',
                      C->size[0], C->size[1], A.size[1],
                      alpha,
                      A.data, A.size[0],
                      B.data, B.size[0],
                      beta,
                      C->data, C->size[0]);
    break;
  case 1:
    if (A.size[0] != B.size[0]) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_matmul",
                             "A and B have incompatible dimensions");
    }
    C->size[0] = A.size[1];
    C->size[1] = B.size[1];
    rc = qmckl_dgemm (context, 'T', 'N',
                      C->size[0], C->size[1], A.size[0],
                      alpha,
                      A.data, A.size[0],
                      B.data, B.size[0],
                      beta,
                      C->data, C->size[0]);
    break;
  case 2:
    if (A.size[1] != B.size[1]) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_matmul",
                             "A and B have incompatible dimensions");
    }
    C->size[0] = A.size[0];
    C->size[1] = B.size[0];
    rc = qmckl_dgemm (context, 'N', 'T',
                      C->size[0], C->size[1], A.size[1],
                      alpha,
                      A.data, A.size[0],
                      B.data, B.size[0],
                      beta,
                      C->data, C->size[0]);
    break;
  case 3:
    if (A.size[0] != B.size[1]) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_matmul",
                             "A and B have incompatible dimensions");
    }
    C->size[0] = A.size[1];
    C->size[1] = B.size[0];
    rc = qmckl_dgemm (context, 'T', 'T',
                      C->size[0], C->size[1], A.size[0],
                      alpha,
                      A.data, A.size[0],
                      B.data, B.size[0],
                      beta,
                      C->data, C->size[0]);
    break;
  }
  return rc;
}

qmckl_exit_code
qmckl_transpose (qmckl_context context,
                 const qmckl_matrix A,
                 qmckl_matrix At )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (A.size[0] < 1) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_transpose",
                           "Invalid size for A");
  }

  if (At.data == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_transpose",
                           "Output matrix not allocated");
  }

  if (At.size[0] != A.size[1] || At.size[1] != A.size[0]) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_transpose",
                           "Invalid size for At");
  }

  for (int64_t j=0 ; j<At.size[1] ; ++j)
    for (int64_t i=0 ; i<At.size[0] ; ++i)
      qmckl_mat(At, i, j) = qmckl_mat(A, j, i);

  return QMCKL_SUCCESS;
}

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_transpose_device (qmckl_context context,
                        const qmckl_matrix A,
                        qmckl_matrix At )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (A.size[0] < 1) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_transpose",
                           "Invalid size for A");
  }

  if (At.data == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_transpose",
                           "Output matrix not allocated");
  }

  if (At.size[0] != A.size[1] || At.size[1] != A.size[0]) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_transpose",
                           "Invalid size for At");
  }

  double * A_data = A.data_device;
  int A_s0 = A.size[0];

  double * At_data = At.data_device;
  int At_s0 = At.size[0];
  int At_s1 = At.size[1];

  #pragma omp target is_device_ptr(A_data, At_data)
  {
  #pragma omp parallel for collapse(2)
  for (int64_t j=0 ; j<At_s1 ; ++j)
    for (int64_t i=0 ; i<At_s0 ; ++i)
      At_data[i + j*At_s0] = A_data[j + i*A_s0];
  }

  return QMCKL_SUCCESS;
}
#endif
