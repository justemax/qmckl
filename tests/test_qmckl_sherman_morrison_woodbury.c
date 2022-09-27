/* [[file:~/Documents/work/qmckl_myFork/qmckl/org/qmckl_sherman_morrison_woodbury.org::*Headers][Headers:2]] */
#include "qmckl.h"
#include "assert.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>

int main() {
  qmckl_context context;
  context = qmckl_context_create();
  qmckl_exit_code rc;
/* Headers:2 ends here */

const uint64_t Dim = 21;
const double breakdown = 1e-3;
const double tolerance = 1e-3;
double res[441];

#include "sm_test.h"

assert(Updates1 != NULL);
assert(Updates_index1 != NULL);
assert(Slater_inv1 != NULL);

// original determinant of Slater1 (before applying updates)
double det = 3.407025646103221e-10;
rc = qmckl_sherman_morrison(context, Dim, Dim, N_updates1, Updates1, Updates_index1, breakdown, Slater_inv1, &det);

// Check that the determinant is updated properly
assert(fabs(det + 4.120398385068217e-10) < 1e-15);

for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    res[i * Dim + j] = 0;
    for (unsigned int k = 0; k < Dim; k++) {
      res[i * Dim + j] += Slater1[i * Dim + k] * Slater_inv1[k * Dim + j];
    }
  }
}
rc = QMCKL_SUCCESS;
for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    if (i == j && fabs(res[i * Dim + j] - 1) > tolerance) {
      rc = QMCKL_FAILURE;
    }
    if (i != j && fabs(res[i * Dim + j]) > tolerance) {
      rc = QMCKL_FAILURE;
    }
  }
}
assert(rc == QMCKL_SUCCESS);

assert(Updates2 != NULL);
assert(Updates_index2 != NULL);
assert(Slater_inv2 != NULL);
det = -1.4432116661319376e-11;
rc = qmckl_woodbury_2(context, Dim, Dim, Updates2, Updates_index2, breakdown, Slater_inv2, &det);
assert(fabs(det-2.367058141251457e-10) < 1e-15);
for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    res[i * Dim + j] = 0;
    for (unsigned int k = 0; k < Dim; k++) {
      res[i * Dim + j] += Slater2[i * Dim + k] * Slater_inv2[k * Dim + j];
    }
  }
}
rc = QMCKL_SUCCESS;
for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    if (i == j && fabs(res[i * Dim + j] - 1) > tolerance) {
      rc = QMCKL_FAILURE;
    }
    if (i != j && fabs(res[i * Dim + j]) > tolerance) {
      rc = QMCKL_FAILURE;
    }
  }
}
assert(rc == QMCKL_SUCCESS);

assert(Updates3 != NULL);
assert(Updates_index3 != NULL);
assert(Slater_inv3_1 != NULL);
det = -1.23743195512859e-09;
rc = qmckl_woodbury_3(context, Dim, Dim, Updates3, Updates_index3, breakdown, Slater_inv3_1, &det);
assert(fabs(det - 1.602708950725074e-10) < 1e-15);
for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    res[i * Dim + j] = 0;
    for (unsigned int k = 0; k < Dim; k++) {
      res[i * Dim + j] += Slater3[i * Dim + k] * Slater_inv3_1[k * Dim + j];
    }
  }
}
rc = QMCKL_SUCCESS;
for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    if (i == j && fabs(res[i * Dim + j] - 1) > tolerance) {
      rc = QMCKL_FAILURE;
    }
    if (i != j && fabs(res[i * Dim + j]) > tolerance) {
      rc = QMCKL_FAILURE;
    }
  }
}
assert(rc == QMCKL_SUCCESS);

assert(Updates3 != NULL);
assert(Updates_index3 != NULL);
assert(Slater_inv3_2 != NULL);
det = -1.23743195512859e-09;
rc = qmckl_sherman_morrison_splitting(context, Dim, Dim, N_updates3, Updates3, Updates_index3, breakdown, Slater_inv3_2, &det);
assert(fabs(det - 1.602708950725074e-10) < 1e-15);
for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    res[i * Dim + j] = 0;
    for (unsigned int k = 0; k < Dim; k++) {
      res[i * Dim + j] += Slater3[i * Dim + k] * Slater_inv3_2[k * Dim + j];
    }
  }
}
rc = QMCKL_SUCCESS;
for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    if (i == j && fabs(res[i * Dim + j] - 1) > tolerance) {
      rc = QMCKL_FAILURE;
    }
    if (i != j && fabs(res[i * Dim + j]) > tolerance) {
      rc = QMCKL_FAILURE;
    }
  }
}
assert(rc == QMCKL_SUCCESS);

assert(Updates5 != NULL);
assert(Updates_index5 != NULL);
assert(Slater_inv5 != NULL);
det = -3.186005284713128e-10;
rc = qmckl_sherman_morrison_smw32s(context, Dim, Dim, N_updates5, Updates5, Updates_index5, breakdown, Slater_inv5, &det);
assert(fabs(det + 5.260200118412903e-10) < 1e-15);

for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    res[i * Dim + j] = 0;
    for (unsigned int k = 0; k < Dim; k++) {
      res[i * Dim + j] += Slater5[i * Dim + k] * Slater_inv5[k * Dim + j];
    }
  }
}
rc = QMCKL_SUCCESS;
for (unsigned int i = 0; i < Dim; i++) {
  for (unsigned int j = 0; j < Dim; j++) {
    if (i == j && fabs(res[i * Dim + j] - 1) > tolerance) {
      rc = QMCKL_FAILURE;
    }
    if (i != j && fabs(res[i * Dim + j]) > tolerance) {
      rc = QMCKL_FAILURE;
    }
  }
}
assert(rc == QMCKL_SUCCESS);

/* [[file:~/Documents/work/qmckl_myFork/qmckl/org/qmckl_sherman_morrison_woodbury.org::*End of files][End of files:1]] */
assert (qmckl_context_destroy(context) == QMCKL_SUCCESS);
  return 0;

}
/* End of files:1 ends here */
