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

/* Reference input data */
int64_t walk_num      = chbrclf_walk_num;
int64_t elec_num      = chbrclf_elec_num;
int64_t elec_up_num   = chbrclf_elec_up_num;
int64_t elec_dn_num   = chbrclf_elec_dn_num;
double  rescale_factor_kappa_ee   = 1.0;
double  rescale_factor_kappa_en   = 1.0;
double  nucl_rescale_factor_kappa = 1.0;
double* elec_coord    = &(chbrclf_elec_coord[0][0][0]);

int64_t  nucl_num      = chbrclf_nucl_num;
double*  charge        = chbrclf_charge;
double*  nucl_coord    = &(chbrclf_nucl_coord[0][0]);

/* --- */

qmckl_exit_code rc;

assert(!qmckl_electron_provided(context));

int64_t n;
rc = qmckl_get_electron_num (context, &n);
assert(rc == QMCKL_NOT_PROVIDED);

rc = qmckl_get_electron_up_num (context, &n);
assert(rc == QMCKL_NOT_PROVIDED);

rc = qmckl_get_electron_down_num (context, &n);
assert(rc == QMCKL_NOT_PROVIDED);


rc = qmckl_set_electron_num (context, elec_up_num, elec_dn_num);
assert(rc == QMCKL_SUCCESS);
assert(qmckl_electron_provided(context));

rc = qmckl_get_electron_up_num (context, &n);
assert(rc == QMCKL_SUCCESS);
assert(n == elec_up_num);

rc = qmckl_get_electron_down_num (context, &n);
assert(rc == QMCKL_SUCCESS);
assert(n == elec_dn_num);

rc = qmckl_get_electron_num (context, &n);
assert(rc == QMCKL_SUCCESS);
assert(n == elec_num);

double k_ee = 0.;
double k_en = 0.;
rc = qmckl_get_electron_rescale_factor_ee (context, &k_ee);
assert(rc == QMCKL_SUCCESS);
assert(k_ee == 1.0);

rc = qmckl_get_electron_rescale_factor_en (context, &k_en);
assert(rc == QMCKL_SUCCESS);
assert(k_en == 1.0);

rc = qmckl_set_electron_rescale_factor_en(context, rescale_factor_kappa_en);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_set_electron_rescale_factor_ee(context, rescale_factor_kappa_ee);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_get_electron_rescale_factor_ee (context, &k_ee);
assert(rc == QMCKL_SUCCESS);
assert(k_ee == rescale_factor_kappa_ee);

rc = qmckl_get_electron_rescale_factor_en (context, &k_en);
assert(rc == QMCKL_SUCCESS);
assert(k_en == rescale_factor_kappa_en);


int64_t w = 0;
rc = qmckl_get_electron_walk_num (context, &w);
assert(rc == QMCKL_SUCCESS);
assert(w == 0);


assert(qmckl_electron_provided(context));

rc = qmckl_set_electron_coord (context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_get_electron_walk_num (context, &w);
assert(rc == QMCKL_SUCCESS);
assert(w == walk_num);

double elec_coord2[walk_num*3*elec_num];

rc = qmckl_get_electron_coord (context, 'N', elec_coord2, walk_num*3*elec_num);
assert(rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i<3*elec_num*walk_num ; ++i) {
  printf("%f %f\n",  elec_coord[i],  elec_coord2[i]);
  assert( elec_coord[i] == elec_coord2[i] );
 }

assert(qmckl_electron_provided(context));


double ee_distance[walk_num * elec_num * elec_num];
rc = qmckl_get_electron_ee_distance(context, ee_distance);

// (e1,e2,w)
// (0,0,0) == 0.
assert(ee_distance[0] == 0.);

// (1,0,0) == (0,1,0)
assert(ee_distance[1] == ee_distance[elec_num]);

// value of (1,0,0)
assert(fabs(ee_distance[1]-7.152322512964209) < 1.e-12);

// (0,0,1) == 0.
assert(ee_distance[elec_num*elec_num] == 0.);

// (1,0,1) == (0,1,1)
assert(ee_distance[elec_num*elec_num+1] == ee_distance[elec_num*elec_num+elec_num]);

// value of (1,0,1)
assert(fabs(ee_distance[elec_num*elec_num+1]-6.5517646321055665) < 1.e-12);

assert(qmckl_electron_provided(context));


double ee_distance_rescaled[walk_num * elec_num * elec_num];
rc = qmckl_get_electron_ee_distance_rescaled(context, ee_distance_rescaled);

// (e1,e2,w)
// (0,0,0) == 0.
assert(ee_distance_rescaled[0] == 0.);

// (1,0,0) == (0,1,0)
assert(ee_distance_rescaled[1] == ee_distance_rescaled[elec_num]);

// value of (1,0,0)
assert(fabs(ee_distance_rescaled[1]-0.9992169566605263) < 1.e-12);

// (0,0,1) == 0.
assert(ee_distance_rescaled[elec_num*elec_num] == 0.);

// (1,0,1) == (0,1,1)
assert(ee_distance_rescaled[elec_num*elec_num+1] == ee_distance_rescaled[elec_num*elec_num+elec_num]);

// value of (1,0,1)
assert(fabs(ee_distance_rescaled[elec_num*elec_num+1]-0.9985724058042633) < 1.e-12);

assert(qmckl_electron_provided(context));


double ee_distance_rescaled_deriv_e[4 * walk_num * elec_num * elec_num];
rc = qmckl_get_electron_ee_distance_rescaled_deriv_e(context, ee_distance_rescaled_deriv_e);

// TODO: Get exact values
//// (e1,e2,w)
//// (0,0,0) == 0.
//assert(ee_distance[0] == 0.);
//
//// (1,0,0) == (0,1,0)
//assert(ee_distance[1] == ee_distance[elec_num]);
//
//// value of (1,0,0)
//assert(fabs(ee_distance[1]-7.152322512964209) < 1.e-12);
//
//// (0,0,1) == 0.
//assert(ee_distance[elec_num*elec_num] == 0.);
//
//// (1,0,1) == (0,1,1)
//assert(ee_distance[elec_num*elec_num+1] == ee_distance[elec_num*elec_num+elec_num]);
//
//// value of (1,0,1)
//assert(fabs(ee_distance[elec_num*elec_num+1]-6.5517646321055665) < 1.e-12);

double ee_potential[walk_num];

rc = qmckl_get_electron_ee_potential(context, &(ee_potential[0]));
assert (rc == QMCKL_SUCCESS);

assert(!qmckl_nucleus_provided(context));
assert(qmckl_electron_provided(context));

rc = qmckl_set_nucleus_num (context, nucl_num);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_set_nucleus_charge (context, charge, nucl_num);
assert (rc == QMCKL_SUCCESS);

rc = qmckl_set_nucleus_coord (context, 'T', nucl_coord, 3*nucl_num);
assert (rc == QMCKL_SUCCESS);

assert(qmckl_nucleus_provided(context));

double en_distance[walk_num][nucl_num][elec_num];

rc = qmckl_get_electron_en_distance(context, &(en_distance[0][0][0]));
assert (rc == QMCKL_SUCCESS);

// (e,n,w) in Fortran notation
// (1,1,1)
assert(fabs(en_distance[0][0][0] - 7.546738741619978) < 1.e-12);

// (1,2,1)
assert(fabs(en_distance[0][1][0] - 8.77102435246984) < 1.e-12);

// (2,1,1)
assert(fabs(en_distance[0][0][1] - 3.698922010513608) < 1.e-12);

// (1,1,2)
assert(fabs(en_distance[1][0][0] - 5.824059436060509) < 1.e-12);

// (1,2,2)
assert(fabs(en_distance[1][1][0] - 7.080482110317645) < 1.e-12);

// (2,1,2)
assert(fabs(en_distance[1][0][1] - 3.1804527583077356) < 1.e-12);

assert(qmckl_electron_provided(context));
assert(qmckl_nucleus_provided(context));

double en_distance_rescaled[walk_num][nucl_num][elec_num];

rc = qmckl_get_electron_en_distance_rescaled(context, &(en_distance_rescaled[0][0][0]));

assert (rc == QMCKL_SUCCESS);

// (e,n,w) in Fortran notation
// (1,1,1)
assert(fabs(en_distance_rescaled[0][0][0] - 0.9994721712909764) < 1.e-12);

// (1,2,1)
assert(fabs(en_distance_rescaled[0][1][0] - 0.9998448354439821) < 1.e-12);

// (2,1,1)
assert(fabs(en_distance_rescaled[0][0][1] - 0.9752498074577688) < 1.e-12);

// (1,1,2)
assert(fabs(en_distance_rescaled[1][0][0] - 0.9970444172399963) < 1.e-12);

// (1,2,2)
assert(fabs(en_distance_rescaled[1][1][0] - 0.9991586325813303) < 1.e-12);

// (2,1,2)
assert(fabs(en_distance_rescaled[1][0][1] - 0.9584331688679852) < 1.e-12);

assert(qmckl_electron_provided(context));

rc = qmckl_set_nucleus_rescale_factor (context, nucl_rescale_factor_kappa);
assert(rc == QMCKL_SUCCESS);

assert(qmckl_nucleus_provided(context));

double en_distance_rescaled_deriv_e[walk_num][4][nucl_num][elec_num];

rc = qmckl_get_electron_en_distance_rescaled_deriv_e(context, &(en_distance_rescaled_deriv_e[0][0][0][0]));

assert (rc == QMCKL_SUCCESS);

// TODO: check exact values
//// (e,n,w) in Fortran notation
//// (1,1,1)
//assert(fabs(en_distance_rescaled[0][0][0] - 7.546738741619978) < 1.e-12);
//
//// (1,2,1)
//assert(fabs(en_distance_rescaled[0][1][0] - 8.77102435246984) < 1.e-12);
//
//// (2,1,1)
//assert(fabs(en_distance_rescaled[0][0][1] - 3.698922010513608) < 1.e-12);
//
//// (1,1,2)
//assert(fabs(en_distance_rescaled[1][0][0] - 5.824059436060509) < 1.e-12);
//
//// (1,2,2)
//assert(fabs(en_distance_rescaled[1][1][0] - 7.080482110317645) < 1.e-12);
//
//// (2,1,2)
//assert(fabs(en_distance_rescaled[1][0][1] - 3.1804527583077356) < 1.e-12);

double en_potential[walk_num];

rc = qmckl_get_electron_en_potential(context, &(en_potential[0]));
assert (rc == QMCKL_SUCCESS);

if (qmckl_context_destroy(context) != QMCKL_SUCCESS)
    return QMCKL_FAILURE;
  return 0;
}
