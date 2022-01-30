#ifndef CBF_PARAMETERS_H
#define CBF_PARAMETERS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>

#include "CBF_driver.h"
#define N_PARAMETER 172 // number of parameters
#define isAFib 0 // hard coded to healthy (not afib)

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
    realtype Qli, Qlao, Qlo, Q_ao_bra, Q_ao_tho, Q_ao_abd;
    realtype Qupi, Qupo;
    realtype Qsp1, Qsp2;
    // realtype QkRi, QkLi;
    // realtype Q_kidneys[24];
    realtype Qll1, Qll2;
    realtype Qk1, Qk2;

    realtype radius[45];
    realtype lengths[45];
    realtype tau[45]; // shearstress
    realtype* flows[45];
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

    realtype P_cerebral[54];

    // realtype x0_aut[6], x_aut[6], A_CO2[6], x0_CO2[6], x_CO2[6], C_dqs[6], dC_dqs[6];

    // for reading in patient data
    realtype P_a, PetCO2;

} * UserData;

typedef struct {
    double currentTime;
    double value;
} TimeSeriesData;

UserData CBF_parameters();

#endif