/* ~qmckl_distance_sq~ */
/*    :PROPERTIES: */
/*    :Name:     qmckl_distance_sq */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    ~qmckl_distance_sq~ computes the matrix of the squared distances */
/*    between all pairs of points in two sets, one point within each set: */

/*    \[ */
/*    C_{ij} = \sum_{k=1}^3 (A_{k,i}-B_{k,j})^2 */
/*    \] */

/*    #+NAME: qmckl_distance_sq_args */
/*    | Variable  | Type             | In/Out | Description                                   | */
/*    |-----------+------------------+--------+-----------------------------------------------| */
/*    | ~context~ | ~qmckl_context~  | in     | Global state                                  | */
/*    | ~transa~  | ~char~           | in     | Array ~A~ is ~'N'~: Normal, ~'T'~: Transposed | */
/*    | ~transb~  | ~char~           | in     | Array ~B~ is ~'N'~: Normal, ~'T'~: Transposed | */
/*    | ~m~       | ~int64_t~        | in     | Number of points in the first set             | */
/*    | ~n~       | ~int64_t~        | in     | Number of points in the second set            | */
/*    | ~A~       | ~double[][lda]~  | in     | Array containing the $m \times 3$ matrix $A$  | */
/*    | ~lda~     | ~int64_t~        | in     | Leading dimension of array ~A~                | */
/*    | ~B~       | ~double[][ldb]~  | in     | Array containing the $n \times 3$ matrix $B$  | */
/*    | ~ldb~     | ~int64_t~        | in     | Leading dimension of array ~B~                | */
/*    | ~C~       | ~double[n][ldc]~ | out    | Array containing the $m \times n$ matrix $C$  | */
/*    | ~ldc~     | ~int64_t~        | in     | Leading dimension of array ~C~                | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~m > 0~ */
/*     - ~n > 0~ */
/*     - ~lda >= 3~ if ~transa == 'N'~ */
/*     - ~lda >= m~ if ~transa == 'T'~ */
/*     - ~ldb >= 3~ if ~transb == 'N'~ */
/*     - ~ldb >= n~ if ~transb == 'T'~ */
/*     - ~ldc >= m~ */
/*     - ~A~ is allocated with at least $3 \times m \times 8$ bytes */
/*     - ~B~ is allocated with at least $3 \times n \times 8$ bytes */
/*     - ~C~ is allocated with at least $m \times n \times 8$ bytes */

/*     #+CALL: generate_c_header(table=qmckl_distance_sq_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_distance_sq (
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
      const int64_t ldc );

/* C header */

/*     #+CALL: generate_c_header(table=qmckl_distance_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_distance (
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
      const int64_t ldc );

/* Device */


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
);
#endif

/* C header */

/*     #+CALL: generate_c_header(table=qmckl_distance_rescaled_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_distance_rescaled (
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
      const double rescale_factor_kappa );

/* ~qmckl_distance_rescaled_deriv_e~ */
/*    :PROPERTIES: */
/*    :Name:     qmckl_distance_rescaled_deriv_e */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    ~qmckl_distance_rescaled_deriv_e~ computes the matrix of the gradient and laplacian of the */
/*    rescaled distance with respect to the electron coordinates. The derivative is a rank 3 tensor. */
/*    The first dimension has a dimension of 4 of which the first three coordinates */
/*    contains the gradient vector and the last index is the laplacian. */


/*    \[ */
/*    C_{ij} = \left( 1 - \exp{-\kappa C_{ij}}\right)/\kappa */
/*    \] */

/*    Here the gradient is defined as follows: */

/*    \[ */
/*    \nabla (C_{ij}(\mathbf{r}_{ee})) = \left(\frac{\delta C_{ij}(\mathbf{r}_{ee})}{\delta x},\frac{\delta C_{ij}(\mathbf{r}_{ee})}{\delta y},\frac{\delta C_{ij}(\mathbf{r}_{ee})}{\delta z} \right) */
/*    \] */
/*    and the laplacian is defined as follows: */

/*    \[ */
/*    \triangle (C_{ij}(r_{ee})) = \frac{\delta^2}{\delta x^2} + \frac{\delta^2}{\delta y^2} + \frac{\delta^2}{\delta z^2} */
/*    \] */

/*    Using the above three formulae, the expression for the gradient and laplacian is */
/*    as follows: */

/*    \[ */
/*    \frac{\delta  C_{ij}(\mathbf{r}_{ee})}{\delta x} = \frac{|(x_i - x_j)|}{r_{ij}} (1 - \kappa R_{ij}) */
/*    \] */

/*    \[ */
/*    \frac{\delta  C_{ij}(\mathbf{r}_{ee})}{\delta y} = \frac{|(y_i - y_j)|}{r_{ij}} (1 - \kappa R_{ij}) */
/*    \] */

/*    \[ */
/*    \frac{\delta  C_{ij}(\mathbf{r}_{ee})}{\delta z} = \frac{|(z_i - z_j)|}{r_{ij}} (1 - \kappa R_{ij}) */
/*    \] */

/*    \[ */
/*    \Delta(C_{ij}(r_{ee}) = \left[ \frac{2}{r_{ij}} - \kappa  \right] (1-\kappa R_{ij}) */
/*    \] */

/*    If the input array is normal (~'N'~), the xyz coordinates are in */
/*    the leading dimension: ~[n][3]~ in C and ~(3,n)~ in Fortran. */

/*    #+NAME: qmckl_distance_rescaled_deriv_e_args */
/*    | Variable               | Type                | In/Out | Description                                           | */
/*    |------------------------+---------------------+--------+-------------------------------------------------------| */
/*    | ~context~              | ~qmckl_context~     | in     | Global state                                          | */
/*    | ~transa~               | ~char~              | in     | Array ~A~ is ~'N'~: Normal, ~'T'~: Transposed         | */
/*    | ~transb~               | ~char~              | in     | Array ~B~ is ~'N'~: Normal, ~'T'~: Transposed         | */
/*    | ~m~                    | ~int64_t~           | in     | Number of points in the first set                     | */
/*    | ~n~                    | ~int64_t~           | in     | Number of points in the second set                    | */
/*    | ~A~                    | ~double[][lda]~     | in     | Array containing the $m \times 3$ matrix $A$          | */
/*    | ~lda~                  | ~int64_t~           | in     | Leading dimension of array ~A~                        | */
/*    | ~B~                    | ~double[][ldb]~     | in     | Array containing the $n \times 3$ matrix $B$          | */
/*    | ~ldb~                  | ~int64_t~           | in     | Leading dimension of array ~B~                        | */
/*    | ~C~                    | ~double[4][n][ldc]~ | out    | Array containing the $4 \times m \times n$ matrix $C$ | */
/*    | ~ldc~                  | ~int64_t~           | in     | Leading dimension of array ~C~                        | */
/*    | ~rescale_factor_kappa~ | ~double~            | in     | Factor for calculating rescaled distances derivatives | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~m > 0~ */
/*     - ~n > 0~ */
/*     - ~lda >= 3~ if ~transa == 'N'~ */
/*     - ~lda >= m~ if ~transa == 'T'~ */
/*     - ~ldb >= 3~ if ~transb == 'N'~ */
/*     - ~ldb >= n~ if ~transb == 'T'~ */
/*     - ~ldc >= m~ */
/*     - ~A~ is allocated with at least $3 \times m \times 8$ bytes */
/*     - ~B~ is allocated with at least $3 \times n \times 8$ bytes */
/*     - ~C~ is allocated with at least $4 \times m \times n \times 8$ bytes */
      
/*     #+CALL: generate_c_header(table=qmckl_distance_rescaled_deriv_e_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_distance_rescaled_deriv_e (
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
      const double rescale_factor_kappa );
