#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <cvodes/cvodes_dense.h>

// # include "ranlib.h"
// # include "rnglib.h"

#define Ith(v,i)			NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ    */
#define IJth(A,i,j)		DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ    */
#define RTOL			RCONST(1.0e-3)   /* scalar relative tolerance            */
#define ATOL			RCONST(1.0e-6)   /* scalar absolute tolerance components */
#define MAXSTEPS	500000

#define NEQ  			57               /* number of equations, ODEs */
#define N_PARAMETER 172         	// number of parameters
#define DELTAT 		0.001  	// time step, seconds. You cannot have time step bigger than a few milliseconds as there are fast processes going on.
#define debg			1 		// this is always 1.

// a check to choose if we want to output all or few state variables
// 1 means all, 0 means only important pressures.
#define numBeats 	1000 // this is total number of beats in a simulation run

typedef struct {
/* Model parameters, some of which are time dependent. */
realtype Ts;
realtype Tasys;
realtype Tav;
realtype p_ursino[N_PARAMETER];
realtype p_C[27];
realtype p_R[45];

realtype Edias_la, Esys_la, Edias_ra, Esys_ra;
realtype Edias_lv, Esys_lv, Edias_rv, Esys_rv;
realtype Ela, Era, Elv, Erv;
realtype Qsup, Qab, Qinf, Qrao, Qro, Qpa;
realtype Qli, Qlao,	Qlo, Q_ao_bra, Q_ao_tho, Q_ao_abd;
realtype Qupi, Qupo	;
realtype Qsp1, Qsp2;
// realtype QkRi, QkLi;
// realtype Q_kidneys[24];
realtype Qll1, 	Qll2;
realtype Qk1, Qk2;

realtype radius[45];
realtype lengths[45];
realtype tau[45]; //shearstress
realtype *flows[45];
realtype argv[5];
realtype RR[3];
int CoW;

// Cerebral blood flow
// realtype C_dqs[6];
// realtype C_dqs_old[6];
realtype CBF;
realtype q_ml;
realtype q_al;
realtype q_pl;
realtype q_mr;
realtype q_ar;
realtype q_pr;
realtype q_ICAl, q_ICAr, q_PCoAl, q_PCoAr, q_MCAl, q_MCAr, q_ACAl, q_ACAr, q_PCAl, q_PCAr;
realtype temp;

// BAROREFLEX
realtype SNA, PNA, SNA_buffer[5000], PNA_buffer[500];

realtype P_cerebral[54]

// realtype x0_aut[6], x_aut[6], A_CO2[6], x0_CO2[6], x_CO2[6], C_dqs[6], dC_dqs[6];

} *UserData;

#include "f_ursino.c"
// #include "ranlib.c"
// #include "rnglib.c"

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}
