#include "qmckl.h"
#include <assert.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include "n2.h"

int main() {
  qmckl_context context;
  context = qmckl_context_create();

/* Reference input data */
int64_t walk_num      = n2_walk_num;
int64_t elec_num      = n2_elec_num;
int64_t elec_up_num   = n2_elec_up_num;
int64_t elec_dn_num   = n2_elec_dn_num;
double  rescale_factor_kappa_ee   = 1.0;
double  rescale_factor_kappa_en   = 1.0;
double  nucl_rescale_factor_kappa = 1.0;
double* elec_coord    = &(n2_elec_coord[0][0][0]);

const double*   nucl_charge   = n2_charge;
int64_t  nucl_num      = n2_nucl_num;
double*  nucl_coord    = &(n2_nucl_coord[0][0]);
int64_t size_max;

/* Provide Electron data */

qmckl_exit_code rc;

assert(!qmckl_electron_provided(context));

rc = qmckl_set_electron_num (context, elec_up_num, elec_dn_num);
assert(rc == QMCKL_SUCCESS);

assert(qmckl_electron_provided(context));

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


rc = qmckl_set_electron_coord (context, 'N', walk_num, elec_coord, walk_num*3*elec_num);
assert(rc == QMCKL_SUCCESS);

double elec_coord2[walk_num*3*elec_num];

rc = qmckl_get_electron_coord (context, 'N', elec_coord2, walk_num*3*elec_num);
assert(rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i<3*elec_num ; ++i) {
  assert( elec_coord[i] == elec_coord2[i] );
 }


/* Provide Nucleus data */

assert(!qmckl_nucleus_provided(context));

rc = qmckl_set_nucleus_num (context, nucl_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_nucleus_provided(context));

double k;
rc = qmckl_get_nucleus_rescale_factor (context, &k);
assert(rc == QMCKL_SUCCESS);
assert(k == 1.0);


rc = qmckl_set_nucleus_rescale_factor (context, nucl_rescale_factor_kappa);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_get_nucleus_rescale_factor (context, &k);
assert(rc == QMCKL_SUCCESS);
assert(k == nucl_rescale_factor_kappa);

double nucl_coord2[3*nucl_num];

rc = qmckl_get_nucleus_coord (context, 'T', nucl_coord2, 3*nucl_num);
assert(rc == QMCKL_NOT_PROVIDED);

rc = qmckl_set_nucleus_coord (context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

assert(!qmckl_nucleus_provided(context));

rc = qmckl_get_nucleus_coord (context, 'N', nucl_coord2, nucl_num*3);
assert(rc == QMCKL_SUCCESS);
for (int64_t k=0 ; k<3 ; ++k) {
  for (int64_t i=0 ; i<nucl_num ; ++i) {
    assert( nucl_coord[nucl_num*k+i] == nucl_coord2[3*i+k] );
  }
}

rc = qmckl_get_nucleus_coord (context, 'T', nucl_coord2, nucl_num*3);
assert(rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i<3*nucl_num ; ++i) {
  assert( nucl_coord[i] == nucl_coord2[i] );
}

double nucl_charge2[nucl_num];

rc = qmckl_get_nucleus_charge(context, nucl_charge2, nucl_num);
assert(rc == QMCKL_NOT_PROVIDED);

rc = qmckl_set_nucleus_charge(context, nucl_charge, nucl_num);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_get_nucleus_charge(context, nucl_charge2, nucl_num);
assert(rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i<nucl_num ; ++i) {
  assert( nucl_charge[i] == nucl_charge2[i] );
 }
assert(qmckl_nucleus_provided(context));

assert(qmckl_electron_provided(context));

int64_t type_nucl_num = n2_type_nucl_num;
int64_t* type_nucl_vector = &(n2_type_nucl_vector[0]);
int64_t aord_num = n2_aord_num;
int64_t bord_num = n2_bord_num;
int64_t cord_num = n2_cord_num;
double* aord_vector = &(n2_aord_vector[0][0]);
double* bord_vector = &(n2_bord_vector[0]);
double* cord_vector = &(n2_cord_vector[0][0]);
int64_t dim_cord_vect=0;

/* Initialize the Jastrow data */
rc = qmckl_init_jastrow(context);
assert(!qmckl_jastrow_provided(context));

/* Set the data */
rc = qmckl_set_jastrow_ord_num(context, aord_num, bord_num, cord_num);
assert(rc == QMCKL_SUCCESS);
rc = qmckl_set_jastrow_type_nucl_num(context, type_nucl_num);
assert(rc == QMCKL_SUCCESS);
rc = qmckl_set_jastrow_type_nucl_vector(context, type_nucl_vector, nucl_num);
assert(rc == QMCKL_SUCCESS);
rc = qmckl_set_jastrow_aord_vector(context, aord_vector,(aord_num+1)*type_nucl_num);
assert(rc == QMCKL_SUCCESS);
rc = qmckl_set_jastrow_bord_vector(context, bord_vector,(bord_num+1));
assert(rc == QMCKL_SUCCESS);
rc = qmckl_get_jastrow_dim_cord_vect(context, &dim_cord_vect);
assert(rc == QMCKL_SUCCESS);
rc = qmckl_set_jastrow_cord_vector(context, cord_vector,dim_cord_vect*type_nucl_num);
assert(rc == QMCKL_SUCCESS);

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_provided(context));

double asymp_jasb[2];
rc = qmckl_get_jastrow_asymp_jasb(context, asymp_jasb,2);

// calculate asymp_jasb
assert(fabs(asymp_jasb[0]-0.5323750557252571) < 1.e-12);
assert(fabs(asymp_jasb[1]-0.31567342786262853) < 1.e-12);

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_provided(context));

double factor_ee[walk_num];
rc = qmckl_get_jastrow_factor_ee(context, factor_ee, walk_num);

// calculate factor_ee
assert(fabs(factor_ee[0]+4.282760865958113) < 1.e-12);

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_provided(context));

// calculate factor_ee_deriv_e
double factor_ee_deriv_e[walk_num][4][elec_num];
rc = qmckl_get_jastrow_factor_ee_deriv_e(context, &(factor_ee_deriv_e[0][0][0]),walk_num*4*elec_num);

// check factor_ee_deriv_e
assert(fabs(factor_ee_deriv_e[0][0][0]-0.16364894652107934) < 1.e-12);
assert(fabs(factor_ee_deriv_e[0][1][0]+0.6927548119830084 ) < 1.e-12);
assert(fabs(factor_ee_deriv_e[0][2][0]-0.073267755223968  ) < 1.e-12);
assert(fabs(factor_ee_deriv_e[0][3][0]-1.5111672803213185 ) < 1.e-12);

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_provided(context));

double factor_en[walk_num];
rc = qmckl_get_jastrow_factor_en(context, factor_en,walk_num);

// calculate factor_en
assert(fabs(factor_en[0]+5.865822569188727) < 1.e-12);

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_provided(context));

// calculate factor_en_deriv_e
double factor_en_deriv_e[walk_num][4][elec_num];
rc = qmckl_get_jastrow_factor_en_deriv_e(context, &(factor_en_deriv_e[0][0][0]),walk_num*4*elec_num);

// check factor_en_deriv_e
assert(fabs(factor_en_deriv_e[0][0][0]-0.11609919541763383) < 1.e-12);
assert(fabs(factor_en_deriv_e[0][1][0]+0.23301394780804574) < 1.e-12);
assert(fabs(factor_en_deriv_e[0][2][0]-0.17548337641865783) < 1.e-12);
assert(fabs(factor_en_deriv_e[0][3][0]+0.9667363412285741 ) < 1.e-12);

assert(qmckl_electron_provided(context));


double een_rescaled_e[walk_num][(cord_num + 1)][elec_num][elec_num];
rc = qmckl_get_jastrow_een_rescaled_e(context, &(een_rescaled_e[0][0][0][0]),elec_num*elec_num*(cord_num+1)*walk_num);

// value of (0,2,1)
assert(fabs(een_rescaled_e[0][1][0][2]-0.08084493981483197)   < 1.e-12);
assert(fabs(een_rescaled_e[0][1][0][3]-0.1066745707571846)    < 1.e-12);
assert(fabs(een_rescaled_e[0][1][0][4]-0.01754273169464735)   < 1.e-12);
assert(fabs(een_rescaled_e[0][2][1][3]-0.02214680362033448)   < 1.e-12);
assert(fabs(een_rescaled_e[0][2][1][4]-0.0005700154999202759) < 1.e-12);
assert(fabs(een_rescaled_e[0][2][1][5]-0.3424402276009091)    < 1.e-12);

//assert(qmckl_electron_provided(context));
double een_rescaled_e_deriv_e[walk_num][(cord_num + 1)][elec_num][4][elec_num];
size_max=walk_num*(cord_num + 1)*elec_num*4*elec_num;
rc = qmckl_get_jastrow_een_rescaled_e_deriv_e(context,
          &(een_rescaled_e_deriv_e[0][0][0][0][0]),size_max);

// value of (0,0,0,2,1)
assert(fabs(een_rescaled_e_deriv_e[0][1][0][0][2] + 0.05991352796887283   ) < 1.e-12);
assert(fabs(een_rescaled_e_deriv_e[0][1][0][0][3] + 0.011714035071545248  ) < 1.e-12);
assert(fabs(een_rescaled_e_deriv_e[0][1][0][0][4] + 0.00441398875758468   ) < 1.e-12);
assert(fabs(een_rescaled_e_deriv_e[0][2][1][0][3] + 0.013553180060167595  ) < 1.e-12);
assert(fabs(een_rescaled_e_deriv_e[0][2][1][0][4] + 0.00041342909359870457) < 1.e-12);
assert(fabs(een_rescaled_e_deriv_e[0][2][1][0][5] + 0.5880599146214673    ) < 1.e-12);

assert(qmckl_electron_provided(context));

double een_rescaled_n[walk_num][(cord_num + 1)][nucl_num][elec_num];
size_max=walk_num*(cord_num + 1)*nucl_num*elec_num;
rc = qmckl_get_jastrow_een_rescaled_n(context, &(een_rescaled_n[0][0][0][0]),size_max);

// value of (0,2,1)
assert(fabs(een_rescaled_n[0][1][0][2]-0.10612983920006765)  < 1.e-12);
assert(fabs(een_rescaled_n[0][1][0][3]-0.135652809635553)    < 1.e-12);
assert(fabs(een_rescaled_n[0][1][0][4]-0.023391817607642338) < 1.e-12);
assert(fabs(een_rescaled_n[0][2][1][3]-0.880957224822116)    < 1.e-12);
assert(fabs(een_rescaled_n[0][2][1][4]-0.027185942659395074) < 1.e-12);
assert(fabs(een_rescaled_n[0][2][1][5]-0.01343938025140174)  < 1.e-12);

assert(qmckl_electron_provided(context));

double een_rescaled_n_deriv_e[walk_num][(cord_num + 1)][nucl_num][4][elec_num];
size_max=walk_num*(cord_num + 1)*nucl_num*4*elec_num;
rc = qmckl_get_jastrow_een_rescaled_n_deriv_e(context, &(een_rescaled_n_deriv_e[0][0][0][0][0]),size_max);

// value of (0,2,1)
assert(fabs(een_rescaled_n_deriv_e[0][1][0][0][2]+0.07633444246999128   )  < 1.e-12);
assert(fabs(een_rescaled_n_deriv_e[0][1][0][0][3]-0.00033282346259738276)  < 1.e-12);
assert(fabs(een_rescaled_n_deriv_e[0][1][0][0][4]+0.004775370547333061  )  < 1.e-12);
assert(fabs(een_rescaled_n_deriv_e[0][2][1][0][3]-0.1362654644223866    )  < 1.e-12);
assert(fabs(een_rescaled_n_deriv_e[0][2][1][0][4]+0.0231253431662794    )  < 1.e-12);
assert(fabs(een_rescaled_n_deriv_e[0][2][1][0][5]-0.001593334817691633  )  < 1.e-12);

assert(qmckl_electron_provided(context));

double tmp_c[walk_num][cord_num][cord_num+1][nucl_num][elec_num];
rc = qmckl_get_jastrow_tmp_c(context, &(tmp_c[0][0][0][0][0]));

double dtmp_c[walk_num][cord_num][cord_num+1][nucl_num][4][elec_num];
rc = qmckl_get_jastrow_dtmp_c(context, &(dtmp_c[0][0][0][0][0][0]));

printf("%e\n%e\n", tmp_c[0][0][1][0][0], 2.7083473948352403);
assert(fabs(tmp_c[0][0][1][0][0] - 2.7083473948352403) < 1e-12);

printf("%e\n%e\n", tmp_c[0][1][0][0][0],0.237440520852232);
assert(fabs(dtmp_c[0][1][0][0][0][0] - 0.237440520852232) < 1e-12);

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_provided(context));

double factor_een[walk_num];
rc = qmckl_get_jastrow_factor_een(context, &(factor_een[0]),walk_num);

assert(fabs(factor_een[0] + 0.37407972141304213) < 1e-12);

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_provided(context));

double factor_een_deriv_e[4][walk_num][elec_num];
rc = qmckl_get_jastrow_factor_een_deriv_e(context, &(factor_een_deriv_e[0][0][0]),4*walk_num*elec_num);

assert(fabs(factor_een_deriv_e[0][0][0] + 0.0005481671107226865) < 1e-12);

rc = qmckl_context_destroy(context);
    assert (rc == QMCKL_SUCCESS);

    return 0;
}
