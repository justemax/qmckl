#include "qmckl.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "chbrclf.h"

int main() {
  qmckl_context context;
  context = qmckl_context_create();

const double*   nucl_charge   = chbrclf_charge;
const double*   nucl_coord    = &(chbrclf_nucl_coord[0][0]);
const double    nucl_rescale_factor_kappa = 2.0;

/* --- */

qmckl_exit_code rc;

assert(!qmckl_nucleus_provided(context));

int64_t n;
rc = qmckl_get_nucleus_num (context, &n);
assert(rc == QMCKL_NOT_PROVIDED);


rc = qmckl_set_nucleus_num (context, chbrclf_nucl_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_nucleus_provided(context));

rc = qmckl_get_nucleus_num (context, &n);
assert(rc == QMCKL_SUCCESS);
assert(n == chbrclf_nucl_num);

double k;
rc = qmckl_get_nucleus_rescale_factor (context, &k);
assert(rc == QMCKL_SUCCESS);
assert(k == 1.0);


rc = qmckl_set_nucleus_rescale_factor (context, nucl_rescale_factor_kappa);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_get_nucleus_rescale_factor (context, &k);
assert(rc == QMCKL_SUCCESS);
assert(k == nucl_rescale_factor_kappa);

double nucl_coord2[3*chbrclf_nucl_num];

rc = qmckl_get_nucleus_coord (context, 'T', nucl_coord2, 3*chbrclf_nucl_num);
assert(rc == QMCKL_NOT_PROVIDED);

rc = qmckl_set_nucleus_coord (context, 'T', &(nucl_coord[0]), 3*chbrclf_nucl_num);
assert(rc == QMCKL_SUCCESS);

assert(!qmckl_nucleus_provided(context));

rc = qmckl_get_nucleus_coord (context, 'N', nucl_coord2, 3*chbrclf_nucl_num);
assert(rc == QMCKL_SUCCESS);
for (size_t k=0 ; k<3 ; ++k) {
  for (int64_t i=0 ; i<chbrclf_nucl_num ; ++i) {
    assert( nucl_coord[chbrclf_nucl_num*k+i] == nucl_coord2[3*i+k] );
  }
}

rc = qmckl_get_nucleus_coord (context, 'T', nucl_coord2, 3*chbrclf_nucl_num);
assert(rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i<3*chbrclf_nucl_num ; ++i) {
  assert( nucl_coord[i] == nucl_coord2[i] );
}

double nucl_charge2[chbrclf_nucl_num];

rc = qmckl_get_nucleus_charge(context, nucl_charge2, chbrclf_nucl_num);
assert(rc == QMCKL_NOT_PROVIDED);

rc = qmckl_set_nucleus_charge(context, nucl_charge, chbrclf_nucl_num);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_get_nucleus_charge(context, nucl_charge2, chbrclf_nucl_num);
assert(rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i<chbrclf_nucl_num ; ++i) {
  assert( nucl_charge[i] == nucl_charge2[i] );
 }
assert(qmckl_nucleus_provided(context));

/* Reference input data */

assert(qmckl_nucleus_provided(context));

double distance[chbrclf_nucl_num*chbrclf_nucl_num];
rc = qmckl_get_nucleus_nn_distance(context, distance, chbrclf_nucl_num*chbrclf_nucl_num);
assert(distance[0] == 0.);
assert(distance[1] == distance[chbrclf_nucl_num]);
assert(fabs(distance[1]-2.070304721365169) < 1.e-12);

/* Reference input data */
/* TODO */

//assert(qmckl_nucleus_provided(context));
//
//double distance[nucl_num*nucl_num];
//rc = qmckl_get_nucleus_nn_distance(context, distance, nucl_num*nucl_num);
//assert(distance[0] == 0.);
//assert(distance[1] == distance[nucl_num]);
//assert(fabs(distance[1]-2.070304721365169) < 1.e-12);

/* Reference input data */

assert(qmckl_nucleus_provided(context));

double rep;
rc = qmckl_get_nucleus_repulsion(context, &rep);
assert(rep - 318.2309879436158 < 1.e-10);

if (qmckl_context_destroy(context) != QMCKL_SUCCESS)
    return QMCKL_FAILURE;
  return 0;
}
