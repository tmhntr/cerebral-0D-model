/*
Timothy J. Hunter, Part of PM3 platforms.
Human circulation, baroreceptor mechanism, and cerebral circulation RHS. Sept 24, 2020.

Sept 24. 2020.
To extend this model using the baro-reflex model.
*/
static int f_ursino(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
    // Assign Y and dY
    realtype Y[NEQ];
    realtype dY[NEQ];

    for (int i = 0; i < NEQ; i++) {
        Y[i] = Ith(y, i + 1);
        dY[i] = 0.0;
    };

    // Create data struct for parameter assignment
    UserData data;
    data = (UserData)user_data;

    /*****************************************************************************************************************************/
    /*****************************************************************************************************************************/
    // %%% ALL PASSED PARAMETERS ARE LISTED BELOW %%%

    double time = data->p_ursino[1];
    double respRate = data->p_ursino[2]; // % respiratory rate; units: (breaths/s)

    double loc_t = data->p_ursino[1] - data->p_ursino[4]; // gnu fmod function. p11 is T_n_1, see above for definition of T_n_1.
    data->p_ursino[3] = loc_t; // passed back to main to calculate cardiac output.

    // State variables are asigned below

    double P_a = Ith(y, 9 + 1);

    /**********************************************************
    ******************* baroreflex ****************************
    **********************************************************/

    double Pbco2 = data->p_ursino[5]; // this is the blood CO2 concentration
    double Pbo2 = data->p_ursino[6]; // this is the blood O2 concentration

    //  Table 3 - Baroreflex control parameters
    //  Afferent compartment:
    double tau_aff = data->p_ursino[7]; //  units: seconds
    double G_aff = data->p_ursino[8];

    //  Central Compartment:
    double tau_c = data->p_ursino[9]; //  Same as tau_aff,  units: s

    //  Efferent compartment:
    double S_p = data->p_ursino[10];
    double PNA_max = data->p_ursino[11]; //  units: Hz
    double PNA_min = data->p_ursino[12]; //  units: Hz
    double S_s = data->p_ursino[13];
    double SNA_max = data->p_ursino[14]; //  units: Hz
    double SNA_min = data->p_ursino[15]; //  units: Hz

    //  Table A3
    double k1 = data->p_ursino[16];
    double k2 = data->p_ursino[17];
    double k3 = data->p_ursino[18];
    double k4 = data->p_ursino[19];
    double k5 = data->p_ursino[20];

    // P_aff
    dY[49] = (-Y[49] + G_aff * (P_a)) / tau_aff;

    //  Central Compartment
    double deltaMAP = 0.0;
    if ((Pbco2 > 40.0) && (Pbo2 < 104.0))
        deltaMAP = k1 + k2 * Pbco2 + k3 / Pbo2;
    else if ((Pbco2 <= 40.0) && (Pbo2 < 104.0))
        deltaMAP = k1 + k2 * 40 + k3 / Pbo2;
    else // if ((Pbco2 > 40.0) && (Pbo2 >= 104.0))
        deltaMAP = k1 + k2 * Pbco2 + k3 / 104.0;

    double P_demand = 90.0 + deltaMAP / 100.0; //  Desrcibed in fig. 2 and on page 793

    // for P_error
    dY[50] = (-Y[50] + 1.0 * (P_demand)) / tau_c; // second term in the P_error eq.

    double P_error = -Y[49] + Y[50]; //  eq. 13 // NOTE Dec 18 changed sign of this term

    //  Efferent Compartment:
    //  eq. 31 chemoreceptor operating point
    double deltaG_SNA;
    if (Pbco2 > 40.0)
        deltaG_SNA = k4 * Pbco2 + k5;
    else
        deltaG_SNA = k4 * 40.0 + k5;

    double k_s = (SNA_max - SNA_min) / (4.0 * S_s); //  eq. 15
    data->SNA = ((SNA_max + SNA_min * exp(P_error / k_s)) / (1.0 + exp(P_error / k_s))) * (1.0 + deltaG_SNA / 100.0); //  eq. 14

    double S_v = S_p; //  paper uses inconsistent notation
    double k_v = (PNA_max - PNA_min) / (4.0 * S_v); //  eq. 17
    data->PNA = (PNA_max - PNA_min * exp(P_error / k_v)) / (1.0 + exp(P_error / k_v)); //  eq. 16

    // Table 3 baroreflex control parameters
    double G_k_s0 = data->p_ursino[21]; // units: beats/min/Hz
    double k_k_s0 = data->p_ursino[22];
    double T_s = data->p_ursino[23]; // units: s
    double G_tau_s0 = 10.0; // beats/min/Hz
    double k_tau_s0 = 0.2;

    double G_v0 = data->p_ursino[24]; // units: beats/min/Hz 45
    double k_v0 = data->p_ursino[25];
    double T_v = data->p_ursino[26]; // units: s

    // Myocardial constriction
    double tau_sigma_lv = data->p_ursino[27]; // units: s
    double T_e_lv = data->p_ursino[28]; // units: s
    double G_eff_lv = data->p_ursino[29]; // units: mmHg/ml/Hz

    double tau_sigma_rv = data->p_ursino[30]; // units: s
    double T_e_rv = data->p_ursino[31]; // units: s
    double G_eff_rv = data->p_ursino[32]; // units: mmHg/ml/Hz

    // Veins
    double tau_sigma_V = data->p_ursino[33]; // units: s
    double T_e_V = data->p_ursino[34]; // units: s
    double G_eff_V = data->p_ursino[35]; // units: ml/Hz

    // Arterioles
    double tau_sigma_R = data->p_ursino[36]; // units: s
    double T_e_R = data->p_ursino[37]; // units: s
    double G_eff_R = data->p_ursino[38]; // units: mmHg/(ml/s)/Hz // EDIT G_eff_R increased from 0.2 to 0.21

    double tau_s;
    if (data->SNA > data->SNA_buffer[0]) {
        double k_tau_s = G_tau_s0 * k_tau_s0;
        tau_s = G_tau_s0 * exp(-k_tau_s * data->SNA); // this is a substituted value based on G_tau_s0 as given eqation goes to zero and crashes
    } else {
        tau_s = data->p_ursino[39];
    }

    double k_s2 = G_k_s0 / k_k_s0; // eq. 22
    double G_s = 1.0 * (1.0 - exp(-k_s2 * data->SNA)); // eq. 21 EDIT gain value changed to 1

    // deltaHR_s
    dY[51] = (-Y[51] + G_s * (data->SNA_buffer[(int)(T_s / DELTAT) - 1])) / tau_s;

    double tau_v;
    if (data->PNA > data->PNA_buffer[0]) {
        tau_v = 0.1;
    } else {
        tau_v = data->p_ursino[40];
    }

    double k_v2 = G_v0 / k_v0; // eq. 27
    double G_v = 1.0 * (1.0 - exp(-k_v2 * data->PNA)); // eq. 26, EDIT gain value changed to 1 for more accurate description of HR used in

    // deltaHR_v
    dY[52] = (-Y[52] + G_v * (data->PNA_buffer[(int)(T_v / DELTAT) - 1])) / tau_v;

    // // SNA Effector Sites
    // again the paper uses inconsistent notation
    double T_sigma_lv = T_e_lv;
    double T_sigma_rv = T_e_rv;
    double T_sigma_V = T_e_V;
    double T_sigma_R = T_e_R;

    // The following are all described by eq. 29
    // Left ventricular contractility
    // sigma_lv
    dY[53] = (-Y[53] + G_eff_lv * (data->SNA_buffer[(int)(T_sigma_lv / DELTAT) - 1])) / tau_sigma_lv;

    // Right ventricular contractility
    // sigma_rv
    dY[54] = (-Y[54] + G_eff_rv * (data->SNA_buffer[(int)(T_sigma_rv / DELTAT) - 1])) / tau_sigma_rv;

    // Venous tone
    // sigma_V
    dY[55] = (-Y[55] + G_eff_V * (data->SNA_buffer[(int)(T_sigma_V / DELTAT) - 1])) / tau_sigma_V;

    // Arterial resistance
    // sigma_R
    dY[56] = (-Y[56] + G_eff_R * (data->SNA_buffer[(int)(T_sigma_R / DELTAT) - 1])) / tau_sigma_R;

    double sigma_lv = Y[53]; /// 1.5; //data->p_ursino[113]; // % units: mmHg/ml
    double sigma_rv = Y[54]; /// 0.935; //data->p_ursino[114]; //  units: mmHg/ml
    // double delta_sigma_V 				= data->p_ursino[115]; //  units: ml
    // double V_gain 							= 0.0; //dY[89]/2785.0; //delta_sigma_V / 2785; //  This is the change in venous volume / total venous volume
    double sigma_R = Y[56]; // data->p_ursino[116]; //  units: mmHg/(ml/s)

    /**********************************************************
    *********************** Heldt *****************************
    **********************************************************/

    // % from Lim 2008, unless it is from Heldt 2002. Notation as close to the two papers as possible.
    // % Flows into and out of each vascular compartment. Find units.
    // % This structure works as a valve; allowing positive flow but not negative.

    /************************ VENA CAVA ************************************/
    // Y[18] is P_ra , Y[10] is P_sup
    if (Y[10] > Y[18]) {
        data->Qsup = (Y[10] - Y[18]) / data->p_R[1];
    } else {
        data->Qsup = 0;
    }

    // Y[4] is P_ab , Y[11] is P_inf
    data->Qab = (Y[4] - Y[11]) / data->p_R[2];

    // Y[11] is P_inf , Y[18] is P_ra
    if (Y[11] > Y[18]) {
        data->Qinf = (Y[11] - Y[18]) / data->p_R[3];
    } else {
        data->Qinf = 0;
    }

    /************************ R.A. OULET into R.V ******************************
    Qrao is the flow coming out of the right atrium,
    Rav is the atrio-ventricular value resistance. */
    // Y[18] is P_ra , Y[12] is P_rv
    if (Y[18] > Y[12]) {
        data->Qrao = (Y[18] - Y[12]) / data->p_R[4];
    } else {
        data->Qrao = 0;
    }

    /************************ R.V. OULET ******************************
    Qro is the flow coming out of the right ventricle,
    Rro is the pulmonary valve resistance.*/
    // Y[12] is P_rv , Y[12] is P_pa
    if (Y[12] > Y[13]) {
        data->Qro = (Y[12] - Y[13]) / data->p_R[5];
    } else {
        data->Qro = 0;
    }

    /********************  PULMONARY ARTERY ***************************************
    Qpa is the flow through the pulmonary artery,
    data->p_R[6] is the pulmonary artery resistance.*/
    // Y[13] is P_pa , Y[14] is P_pv
    data->Qpa = (Y[13] - Y[14]) / data->p_R[6];

    /********************* L.A. INLET / PULMONARY VEIN*****************************
    Qli is the pulmonary venous flow AKA the flow into the L. Atrium,
    data->p_R[7] is the pulmonary vein resistance. */
    // Y[14] is P_pv , Y[17] is P_ra
    if (Y[14] > Y[17]) {
        data->Qli = (Y[14] - Y[17]) / data->p_R[7];
    } else {
        data->Qli = 0;
    }

    /************************ L.A. OULET -> L.V ***********************************
    Qlao is the flow coming out of the left atrium,
    Rmv is the mitral value resistance. */
    // Y[17] is P_ra , Y[8] is P_rv
    if (Y[17] > Y[8]) {
        data->Qlao = (Y[17] - Y[8]) / data->p_R[8];
    } else {
        data->Qlao = 0.0;
    }

    /************************ L.V. OULET / AORTIC *********************************
    Qlo is the flow out of the left ventricle,
    Rlo is the aortic valve resistance. */
    // Y[8] is P_rv , Y[9] is P_a
    if (Y[8] > P_a) {
        data->Qlo = (Y[8] - P_a) / data->p_R[9];
    } else {
        data->Qlo = 0;
    }
    // TODO convert back to one aorta segment
    // data->Q_ao_bra = (Y[12] - Y[50])/data->p_R[10];
    //
    // data->Q_ao_tho = (Y[12] - Y[51])/data->p_R[11];
    //
    // data->Q_ao_abd = (Y[51] - Y[52])/data->p_R[12];

    /*************************** UPPER BODY ***********************************
    Qupi is the inlet flow through the upper body,
    Qupo is the outlet flow through the upper body,
    Rup1 & Rup2, are the resistances of each. */
    // P_a is P_a , Y[0] is P_up
    data->Qupi = (P_a - Y[0]) / (data->p_R[10] * sigma_R);

    // Y[0] is P_up, Y[10] is P_sup
    if (Y[0] > Y[10]) {
        data->Qupo = (Y[0] - Y[10]) / data->p_R[11];
    } else {
        data->Qupo = 0;
    }

    /*************************** SPLANCHNIC  ***********************************
    Qsp1 is the splanchnic microcirculation?
    Qsp2 is the splanchnic venous circulation,
    Rsp1 & Rsp2 are the resistances of each.
    Y[9] is P_a,
    Y[2] is P_sp,
    Y[4] is P_ab
    */
    data->Qsp1 = (P_a - Y[2]) / (data->p_R[12] * sigma_R);

    if (Y[4] < Y[2]) {
        data->Qsp2 = (Y[2] - Y[4]) / data->p_R[13];
    } else {
        data->Qsp2 = 0;
    }

    /*************** DETAILED KIDNEY ************************************/

    // data->QkRi = (Y[12] - Y[48])/data->p_R[17]; // right kidney inlet flow.
    // data->QkLi = (Y[12] - Y[49])/data->p_R[18]; // left kidney inlet flow.
    //
    // data->Qk1 					=  data->QkRi + data->QkLi;
    // data->Qk2           = 0.0;
    // for(i = 0; i<12; i++){
    //   data->Q_kidneys[i]      = (Y[48] - Y[36+i])/(data->p_R[19+i]* sigma_R); // Qk11 to Qk12
    //   data->Q_kidneys[12 + i] = (Y[36+i] - Y[7])/(data->p_R[31+i]* sigma_R);  // Qk21 to Qk212
    //   data->Qk2               = data->Qk2 + data->Q_kidneys[12 + i];
    // }

    /*************** REGULAR KIDNEY ************************************/

    data->Qk1 = (P_a - Y[1]) / data->p_R[14];
    data->Qk2 = (Y[1] - Y[4]) / data->p_R[15];

    /*************************** LOWER BODY ***********************************
    Qll1 & Qll2 are the lower body/legs flow,
    data->p_R[43] & data->p_R[44] are their resistances.
    Y[12] is P_a,
    Y[6] is P_ll,
    Y[7] is P_ab.  */
    data->Qll1 = (P_a - Y[3]) / (data->p_R[16] * sigma_R);

    if (Y[3] > Y[4]) {
        data->Qll2 = (Y[3] - Y[4]) / data->p_R[17];
    } else {
        data->Qll2 = 0;
    }

    /**********************************************************
    *********************** heart ****************************
    **********************************************************/

    double Ts = data->Ts;
    double Tasys = data->Tasys;
    double Tav = data->Tav;
    double act_fxn_v, act_fxn_a;

    /*******************************************************************************
    Ventricular Elastance
    *******************************************************************************/
    if ((loc_t > Tav) && (loc_t <= Tav + Ts)) {
        act_fxn_v = (1.0 - cos(M_PI * (loc_t - Tav) / Ts));
    } else if ((loc_t > Tav + Ts) && (loc_t <= Tav + 1.5 * Ts)) {
        act_fxn_v = (1.0 + cos(2.0 * M_PI * (loc_t - Tav - Ts) / Ts));
    } else {
        act_fxn_v = 0.0;
    }

    double Elv = data->Elv = data->Edias_lv + 0.5 * (data->Esys_lv * sigma_lv - data->Edias_lv) * act_fxn_v;
    double Erv = data->Erv = data->Edias_rv + 0.5 * (data->Esys_rv * sigma_rv - data->Edias_rv) * act_fxn_v;

    double Clv = 1.0 / Elv; //
    double Crv = 1.0 / Erv; //

    /*******************************************************************************
     Atrial Elastance
    *******************************************************************************/
    if ((0 < loc_t) && (loc_t <= Tasys)) {
        act_fxn_a = (1.0 - cos(M_PI * (loc_t / Tasys)));
    } else if ((Tasys <= loc_t) && (loc_t < 1.5 * Tasys)) {
        act_fxn_a = (1.0 + cos(2.0 * M_PI * (loc_t - Tasys) / Tasys));
    } else {
        act_fxn_a = 0.0;
    }

    double Ela = data->Ela = data->Edias_la + 0.5 * (data->Esys_la - data->Edias_la) * act_fxn_a;
    double Era = data->Era = data->Edias_ra + 0.5 * (data->Esys_ra - data->Edias_ra) * act_fxn_a;

    double Cla = 1.0 / Ela;
    double Cra = 1.0 / Era; // new value of elastance ODE variables.
    /**********************************************************
    ********************* cerebral ****************************
    **********************************************************/

    // This cerebral vasculature model is based on the model by Ursino and
    // Giannessi (2010), and includes a detailed circle of willis and cerebral
    // autoregulation and  metabolic reacctivity.
    //
    // Implemented by Timothy Hunter
    // 6 Nov, 2020
    // For PM3 labs
    //
    // Ref: A Model of Cerebrovascular Reactivity Including the Circle
    // of Willis and Cortical Anastomoses, DOI: 10.1007/s10439-010-9923-7

    // This model has 24 ODEs

    // double P_a = Y[9];
    double P_aCO2 = data->p_ursino[5];
    double P_vs = Y[10];

    ///////////////////////// CoW variations //////////////////////////////////////
    // The following are settings that determine which blood vessels are to be excluded from the model
    int A1l, P1l, PCoAl, PCoAr;

    if (data->CoW == 0) {
        // Normal CoW
        A1l = 0;
        P1l = 0;
        PCoAl = 0;
        PCoAr = 0;
    } else if (data->CoW == 1) {
        // Left PCoA missing; PCoA
        A1l = 0;
        P1l = 0;
        PCoAl = 1;
        PCoAr = 0;
    } else if (data->CoW == 2) {
        // Both PCoA missing; PCoAs
        A1l = 0;
        P1l = 0;
        PCoAl = 1;
        PCoAr = 1;
    } else if (data->CoW == 3) {
        // ACAl missing; A1
        A1l = 1;
        P1l = 0;
        PCoAl = 0;
        PCoAr = 0;
    } else if (data->CoW == 4) {
        // PCAl missing; P1
        A1l = 0;
        P1l = 1;
        PCoAl = 0;
        PCoAr = 0;
    } else if (data->CoW == 5) {
        // PCoA and contralateral PCA; PCoA, PCA
        A1l = 0;
        P1l = 1;
        PCoAl = 0;
        PCoAr = 1;
    }

    ///////////////////////////////////////////////////////////////////////////////

    // Parameters from table 1 of ref.
    // Hemodynamic and hydrodynamic
    double R_f = data->P_cerebral[0]; // 2.38 * pow(10,3); // units mmHg s mL^-1
    double R_o = data->P_cerebral[1]; // 526.3; // units mmHg s mL^-1
    double R_pv = data->P_cerebral[2]; // 0.880; // units mmHg s mL^-1
    double R_vs1 = data->P_cerebral[3]; // 0.366; // units mmHg s mL^-1
    double R_cpms = data->P_cerebral[4]; // 120.0; // units mmHg s mL^-1
    double R_cams = data->P_cerebral[5]; // 105.0; // units mmHg s mL^-1
    double R_cpp = data->P_cerebral[6]; // 75.0; // units mmHg s mL^-1
    double R_caa = data->P_cerebral[7]; // 22.0; // units mmHg s mL^-1
    double C_ICAs = data->P_cerebral[8]; // 3.4 * pow(10,-3); // units mL mmHg^-1
    double C_BA = data->P_cerebral[9]; // 1.7 * pow(10,-3); // units mL mmHg^-1
    double r_ICAns = data->P_cerebral[10]; // 0.205; // units cm
    double r_BAn = data->P_cerebral[11]; // 0.17; // units cm
    double r_MCAns = data->P_cerebral[12]; // 0.14; // units cm
    double r_PCA1ns = data->P_cerebral[13]; // 0.1; // units cm
    double r_ACA1ns = data->P_cerebral[14]; // 0.075; // units cm
    double r_PCA2ns = data->P_cerebral[15]; // 0.1; // units cm
    double r_ACA2ns = data->P_cerebral[16]; // 0.075; // units cm
    double r_PCoAns = data->P_cerebral[17]; // 0.036; // units cm
    double r_ACoAn = data->P_cerebral[18]; // 0.04; // units cm
    double l_ICAns = data->P_cerebral[19]; // 13.15; // units cm
    double l_BAn = data->P_cerebral[20]; // 4.92; // units cm
    double l_MCAns = data->P_cerebral[21]; // 7.25; // units cm
    double l_PCA1ns = data->P_cerebral[22]; // 1.0; // units cm
    double l_ACA1ns = data->P_cerebral[23]; // 1.57; // units cm
    double l_PCA2ns = data->P_cerebral[24]; // 4.72; // units cm
    double l_ACA2ns = data->P_cerebral[25]; // 0.672; // units cm
    double l_PCoAns = data->P_cerebral[26]; // 2.0; // units cm
    double l_ACoAn = data->P_cerebral[27]; // 0.5; // units cm
    double k_ven = data->P_cerebral[28]; // 0.155; // units mL^-1
    double k_MCAs = data->P_cerebral[29]; // 12.0;
    double k_E = data->P_cerebral[30]; // 0.077; // units mL^-1
    double P_v1 = data->P_cerebral[31]; // -2.5; // units mmHg
    // P_vs = 6; // units mmHg  NOTE: to be passed as Psup
    double V_dn = data->P_cerebral[32]; // 10.1; // units mL
    double C_dn = data->P_cerebral[33]; // 200.0 * pow(10,-3); // units mL mmHg^-1
    double R_dn = data->P_cerebral[34]; // 5.4; // units mmHg s mL^-1
    double k_R = data->P_cerebral[35]; // 13.1 * pow(10,3); // units mmHg^-3 s mL^-1
    double q_n = data->P_cerebral[36]; // 12.5; // units mL s^-1
    double W_ICAs = data->P_cerebral[37]; // 0.38;
    double W_BA = data->P_cerebral[38]; // 0.24;
    double W_PCA1s = data->P_cerebral[39]; // 0.12;
    double W_ACA1s = data->P_cerebral[40]; // 0.08;
    double W_MCAs = data->P_cerebral[41]; // 0.30;
    double W_PCA2s = data->P_cerebral[42]; // 0.12;
    double W_ACA2s = data->P_cerebral[43]; // 0.08;

    double mu = data->P_cerebral[44]; // 0.004/133.322; // units poise

    // Cerebrlovascular control mechanisms
    double tau_autjs = data->P_cerebral[45]; // 20.0; // units s
    double G_autjs = data->P_cerebral[46]; // 0.9;
    double tau_CO2js = data->P_cerebral[47]; // 40.0; // units s
    double G_CO2js = data->P_cerebral[48]; // 4.0;
    double sat2 = data->P_cerebral[49]; // 7.0;
    double sat1 = data->P_cerebral[50]; // 0.4;
    double k_CO2 = data->P_cerebral[51]; // 15.0;
    double b_CO2 = data->P_cerebral[52]; // 0.5;
    double P_aCO2njs = data->P_cerebral[53]; // 40.0; // units mmHg

    // ****************************************************************************
    //********** compliances ********************

    double C_d0[6], deltaC_d[6], k_Cd[6], C_dqs[6];

    C_d0[0] = C_d0[3] = C_dn * W_MCAs;
    C_d0[1] = C_d0[4] = C_dn * W_ACA2s;
    C_d0[2] = C_d0[5] = C_dn * W_PCA2s;

    for (int i = 0; i < 6; i++) {

        // eq. 25
        if (Y[37 + i] - Y[31 + i] < 0.0) {
            deltaC_d[i] = 2.0 * sat1;
        } else {
            deltaC_d[i] = 2.0 * sat2;
        }

        k_Cd[i] = C_d0[i] * deltaC_d[i] / 4.0;

        C_dqs[i] = C_d0[i] * ((1.0 - deltaC_d[i] / 2.0) + (1.0 + deltaC_d[i] / 2.0) * exp((Y[37 + i] - Y[31 + i]) / k_Cd[i])) / (1.0 + exp((Y[37 + i] - Y[31 + i]) / k_Cd[i]));
        if (C_dqs[i] < 0.0001) {
            C_dqs[i] = 0.0001;
        }
        dY[43 + i] = (C_dqs[i] - Y[43 + i]) / DELTAT;
    }
    //********** end of compliances ********************

    // Derived
    double V_dnms = V_dn * W_MCAs;
    double V_dnas = V_dn * W_ACA2s;
    double V_dnps = V_dn * W_PCA2s;

    double q_nms = q_n * W_MCAs;
    double q_nas = q_n * W_ACA2s;
    double q_nps = q_n * W_PCA2s;

    // C_d0ms = C_dn * W_MCAs;
    // C_d0as = C_dn * W_ACA2s;
    // C_d0ps = C_dn * W_PCA2s;
    //
    // These values are k_R * weight factor
    double k_Rms = k_R / W_MCAs;
    double k_Ras = k_R / W_ACA2s;
    double k_Rps = k_R / W_PCA2s;

    // eq 10
    double R_vs = 0.0;
    if (P_vs < Y[19] & Y[21] != Y[19]) {
        R_vs = (Y[21] - P_vs) * R_vs1 / (Y[21] - Y[19]);
    } else {
        R_vs = R_vs1;
    }

    // eq 2
    double qf = 0.0;
    if (Y[20] > Y[19]) {
        qf = (Y[20] - Y[19]) / R_f;
    } else {
        qf = 0.0;
    }

    // eq 3
    double qo = 0.0;
    if (Y[19] > P_vs) {
        qo = (Y[19] - P_vs) / R_o;
    } else {
        qo = 0.0;
    }

    // eq 9
    double dV_vi = (Y[20] - Y[21]) / R_pv - (Y[21] - P_vs) / R_vs + qo;

    double C_d0ms = C_dn * W_MCAs;
    double C_d0as = C_dn * W_ACA2s;
    double C_d0ps = C_dn * W_PCA2s;

    // eq 26
    double R_dml = (k_Rms * pow(C_d0ms, 2)) / (pow(Y[43] * (Y[22] - Y[19]), 2));
    double R_dal = (k_Ras * pow(C_d0as, 2)) / (pow(Y[44] * (Y[23] - Y[19]), 2));
    double R_dpl = (k_Rps * pow(C_d0ps, 2)) / (pow(Y[45] * (Y[24] - Y[19]), 2));
    double R_dmr = (k_Rms * pow(C_d0ms, 2)) / (pow(Y[46] * (Y[25] - Y[19]), 2));
    double R_dar = (k_Ras * pow(C_d0as, 2)) / (pow(Y[47] * (Y[26] - Y[19]), 2));
    double R_dpr = (k_Rps * pow(C_d0ps, 2)) / (pow(Y[48] * (Y[27] - Y[19]), 2));

    // eq 27
    // These values are either from poiseills formula OR R_dn / Weight factor
    double R_PCA1s = 0.7640; //(8.0 * mu * l_PCA1ns)/(M_PI * pow(r_PCA1ns,4));
    double R_PCA2s = 3.6063; //(92.5 - 58.5) / (q_n * W_PCA2s); //(8.0 * mu * l_PCA2ns)/(M_PI * pow(r_PCA2ns,4)); // R_dn * W_PCA2s;
    double R_PCoAs = 90.9786; //(8.0 * mu * l_PCoAns)/(M_PI * pow(r_PCoAns,4));
    double R_ACA1s = 3.7912; //(8.0 * mu * l_ACA1ns)/(M_PI * pow(r_ACA1ns,4));
    double R_ACA2s = 1.6227; //(92.5 - 58.5) / (q_n * W_ACA2s); // (8.0 * mu * l_ACA2ns)/(M_PI * pow(r_ACA2ns,4));
    double R_ACoA = 14.9228; //(8.0 * mu * l_ACoAn)/(M_PI *pow(r_ACoAn,4));
    double R_ICAs = 0.5689; //(100.0 - 92.5)/(q_n * W_ICAs);   //(8.0 * mu * l_ICAns)/(M_PI * pow(r_ICAns,4)); // R_dn * W_ICAs;
    double R_BA = 0.4501; //(100.0 - 92.5)/(q_n * W_BA);   //(8.0 * mu * l_BAn)/(M_PI * pow(r_BAn,4)); // R_dn * W_BA;

    double P_ICAns = 92.5; // Units: mmHg, ref. p. 962 last paragraph

    // eq 20
    double r_MCAl = r_MCAns * ((1.0 / k_MCAs) * log(Y[28] / P_ICAns) + 1.0);
    double r_MCAr = r_MCAns * ((1.0 / k_MCAs) * log(Y[29] / P_ICAns) + 1.0);

    double R_MCAs = 1.4419; //(8.0 * mu * l_MCAns)/(M_PI * pow(r_MCAl,4));

    // eq 17 (system solved using matlab solve function)
    double P_ACA1l = 0.0;
    double P_ACA1r = 0.0;
    if (A1l == 0) {
        P_ACA1l = (4.0 * Y[28] * R_ACoA * pow(R_ACA2s, 2) + 4.0 * Y[28] * R_ACA1s * pow(R_ACA2s, 2) + 4.0 * Y[29] * R_ACA1s * pow(R_ACA2s, 2) + 4.0 * Y[23] * R_ACoA * pow(R_ACA1s, 2) + 4.0 * Y[23] * pow(R_ACA1s, 2) * R_ACA2s + 4.0 * Y[26] * pow(R_ACA1s, 2) * R_ACA2s + 2.0 * Y[23] * pow(R_ACA1s, 2) * R_dar + 2.0 * Y[26] * pow(R_ACA1s, 2) * R_dal + 4.0 * Y[28] * R_ACoA * R_ACA1s * R_ACA2s + 4.0 * Y[23] * R_ACoA * R_ACA1s * R_ACA2s + 2.0 * Y[28] * R_ACoA * R_ACA1s * R_dal + 2.0 * Y[28] * R_ACoA * R_ACA2s * R_dal + 2.0 * Y[28] * R_ACA1s * R_ACA2s * R_dal + 2.0 * Y[28] * R_ACoA * R_ACA2s * R_dar + 2.0 * Y[28] * R_ACA1s * R_ACA2s * R_dar + 2.0 * Y[29] * R_ACA1s * R_ACA2s * R_dal + 2.0 * Y[29] * R_ACA1s * R_ACA2s * R_dar + 2.0 * Y[23] * R_ACoA * R_ACA1s * R_dar + Y[28] * R_ACoA * R_dal * R_dar + Y[28] * R_ACA1s * R_dal * R_dar + Y[29] * R_ACA1s * R_dal * R_dar) / (4.0 * R_ACoA * pow(R_ACA1s, 2) + 4.0 * R_ACoA * pow(R_ACA2s, 2) + 8.0 * R_ACA1s * pow(R_ACA2s, 2) + 8.0 * pow(R_ACA1s, 2) * R_ACA2s + 2.0 * pow(R_ACA1s, 2) * R_dal + 2.0 * pow(R_ACA1s, 2) * R_dar + 8.0 * R_ACoA * R_ACA1s * R_ACA2s + 2.0 * R_ACoA * R_ACA1s * R_dal + 2.0 * R_ACoA * R_ACA2s * R_dal + 2.0 * R_ACoA * R_ACA1s * R_dar + 4.0 * R_ACA1s * R_ACA2s * R_dal + 2.0 * R_ACoA * R_ACA2s * R_dar + 4.0 * R_ACA1s * R_ACA2s * R_dar + R_ACoA * R_dal * R_dar + 2.0 * R_ACA1s * R_dal * R_dar);

        P_ACA1r = (4.0 * Y[28] * R_ACA1s * pow(R_ACA2s, 2) + 4.0 * Y[29] * R_ACoA * pow(R_ACA2s, 2) + 4.0 * Y[29] * R_ACA1s * pow(R_ACA2s, 2) + 4.0 * Y[23] * pow(R_ACA1s, 2) * R_ACA2s + 4.0 * Y[26] * R_ACoA * pow(R_ACA1s, 2) + 4.0 * Y[26] * pow(R_ACA1s, 2) * R_ACA2s + 2.0 * Y[23] * pow(R_ACA1s, 2) * R_dar + 2.0 * Y[26] * pow(R_ACA1s, 2) * R_dal + 4.0 * Y[29] * R_ACoA * R_ACA1s * R_ACA2s + 4.0 * Y[26] * R_ACoA * R_ACA1s * R_ACA2s + 2.0 * Y[28] * R_ACA1s * R_ACA2s * R_dal + 2.0 * Y[29] * R_ACoA * R_ACA2s * R_dal + 2.0 * Y[28] * R_ACA1s * R_ACA2s * R_dar + 2.0 * Y[29] * R_ACoA * R_ACA1s * R_dar + 2.0 * Y[29] * R_ACA1s * R_ACA2s * R_dal + 2.0 * Y[29] * R_ACoA * R_ACA2s * R_dar + 2.0 * Y[29] * R_ACA1s * R_ACA2s * R_dar + 2.0 * Y[26] * R_ACoA * R_ACA1s * R_dal + Y[28] * R_ACA1s * R_dal * R_dar + Y[29] * R_ACoA * R_dal * R_dar + Y[29] * R_ACA1s * R_dal * R_dar) / (4.0 * R_ACoA * pow(R_ACA1s, 2) + 4.0 * R_ACoA * pow(R_ACA2s, 2) + 8.0 * R_ACA1s * pow(R_ACA2s, 2) + 8.0 * pow(R_ACA1s, 2) * R_ACA2s + 2.0 * pow(R_ACA1s, 2) * R_dal + 2.0 * pow(R_ACA1s, 2) * R_dar + 8.0 * R_ACoA * R_ACA1s * R_ACA2s + 2.0 * R_ACoA * R_ACA1s * R_dal + 2.0 * R_ACoA * R_ACA2s * R_dal + 2.0 * R_ACoA * R_ACA1s * R_dar + 4.0 * R_ACA1s * R_ACA2s * R_dal + 2.0 * R_ACoA * R_ACA2s * R_dar + 4.0 * R_ACA1s * R_ACA2s * R_dar + R_ACoA * R_dal * R_dar + 2.0 * R_ACA1s * R_dal * R_dar);

    } else if (A1l == 1) {
        P_ACA1l = (4.0 * Y[29] * pow(R_ACA2s, 2) + 4.0 * Y[23] * R_ACoA * R_ACA1s + 4.0 * Y[23] * R_ACoA * R_ACA2s + 4.0 * Y[23] * R_ACA1s * R_ACA2s + 4.0 * Y[26] * R_ACA1s * R_ACA2s + 2.0 * Y[29] * R_ACA2s * R_dal + 2.0 * Y[29] * R_ACA2s * R_dar + 2.0 * Y[23] * R_ACoA * R_dar + 2.0 * Y[23] * R_ACA1s * R_dar + 2.0 * Y[26] * R_ACA1s * R_dal + Y[29] * R_dal * R_dar) / (4.0 * R_ACoA * R_ACA1s + 4.0 * R_ACoA * R_ACA2s + 8.0 * R_ACA1s * R_ACA2s + 2.0 * R_ACA1s * R_dal + 2.0 * R_ACoA * R_dar + 2.0 * R_ACA2s * R_dal + 2.0 * R_ACA1s * R_dar + 2.0 * R_ACA2s * R_dar + R_dal * R_dar + 4.0 * pow(R_ACA2s, 2));

        P_ACA1r = (4.0 * Y[29] * pow(R_ACA2s, 2) + 4.0 * Y[29] * R_ACoA * R_ACA2s + 4.0 * Y[23] * R_ACA1s * R_ACA2s + 4.0 * Y[26] * R_ACoA * R_ACA1s + 4.0 * Y[26] * R_ACA1s * R_ACA2s + 2.0 * Y[29] * R_ACoA * R_dar + 2.0 * Y[29] * R_ACA2s * R_dal + 2.0 * Y[29] * R_ACA2s * R_dar + 2.0 * Y[23] * R_ACA1s * R_dar + 2.0 * Y[26] * R_ACA1s * R_dal + Y[29] * R_dal * R_dar) / (4.0 * R_ACoA * R_ACA1s + 4.0 * R_ACoA * R_ACA2s + 8.0 * R_ACA1s * R_ACA2s + 2.0 * R_ACA1s * R_dal + 2.0 * R_ACoA * R_dar + 2.0 * R_ACA2s * R_dal + 2.0 * R_ACA1s * R_dar + 2.0 * R_ACA2s * R_dar + R_dal * R_dar + 4.0 * pow(R_ACA2s, 2));
    }

    // P_ACA1l = (P_dal/(R_ACA2s + R_dal/2) + P_ICAl/R_ACA1s + (P_dar/(R_ACA2s + R_dar/2) + P_ICAr/R_ACA1s) / (1/(R_ACA2s + R_dar/2) + 1/R_ACA1s + 1/R_ACoA)) / ((1/(R_ACA2s + R_dal/2) + 1/R_ACA1s + 1/R_ACoA) - (1/R_ACoA) * (1/(R_ACA2s + R_dar/2) + 1/R_ACA1s + 1/R_ACoA));
    // P_ACA1r = (P_dar/(R_ACA2s + R_dar/2) + P_ICAr/R_ACA1s + (P_dal/(R_ACA2s + R_dal/2) + P_ICAl/R_ACA1s) / (1/(R_ACA2s + R_dal/2) + 1/R_ACA1s + 1/R_ACoA)) / ((1/(R_ACA2s + R_dar/2) + 1/R_ACA1s + 1/R_ACoA) - (1/R_ACoA) * (1/(R_ACA2s + R_dal/2) + 1/R_ACA1s + 1/R_ACoA));
    //
    // double P_ACA1l = (Y[23]/(R_ACA2s + R_dal/2.0) + Y[28]/R_ACA1s + (Y[26]/(R_ACA2s + R_dar/2.0) + Y[29]/R_ACA1s)/(1.0/(R_ACA2s + R_dar/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA)) / (1.0/(R_ACA2s + R_dal/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA - 1.0/(R_ACoA*(1.0/(R_ACA2s + R_dar/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA)));
    // double P_ACA1r = (Y[26]/(R_ACA2s + R_dar/2.0) + Y[29]/R_ACA1s + (Y[23]/(R_ACA2s + R_dal/2.0) + Y[28]/R_ACA1s)/(1.0/(R_ACA2s + R_dal/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA)) / (1.0/(R_ACA2s + R_dar/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA - 1.0/(R_ACoA*(1.0/(R_ACA2s + R_dal/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA)));

    // eq 16
    double P_PCA1l = 0.0;
    if (P1l == 0) {
        if (PCoAl == 0) {
            P_PCA1l = (Y[28] / R_PCoAs + Y[30] / R_PCA1s + Y[24] / (R_PCA2s + R_dpl / 2.0)) / (1.0 / R_PCoAs + 1.0 / R_PCA1s + 1.0 / (R_PCA2s + R_dpl / 2.0));
        } else if (PCoAl == 1) {
            P_PCA1l = (Y[30] / R_PCA1s + Y[24] / (R_PCA2s + R_dpl / 2.0)) / (1.0 / R_PCA1s + 1.0 / (R_PCA2s + R_dpl / 2.0));
        }
    } else if (P1l == 1) {
        if (PCoAl == 0) {
            P_PCA1l = (Y[28] / R_PCoAs + Y[24] / (R_PCA2s + R_dpl / 2.0)) / (1.0 / R_PCoAs + 1.0 / (R_PCA2s + R_dpl / 2.0));
        } else if (PCoAl == 1) {
            P_PCA1l = (Y[24] / (R_PCA2s + R_dpl / 2.0)) / (1.0 / (R_PCA2s + R_dpl / 2.0));
        }
    }
    double P_PCA1r = 0.0;
    if (PCoAr == 0) {
        P_PCA1r = (Y[29] / R_PCoAs + Y[30] / R_PCA1s + Y[27] / (R_PCA2s + R_dpr / 2.0)) / (1.0 / R_PCoAs + 1.0 / R_PCA1s + 1.0 / (R_PCA2s + R_dpr / 2.0));
    } else if (PCoAr == 1) {
        P_PCA1r = (Y[30] / R_PCA1s + Y[27] / (R_PCA2s + R_dpr / 2.0)) / (1.0 / R_PCA1s + 1.0 / (R_PCA2s + R_dpr / 2.0));
    }

    // Flows
    data->q_ICAl = (P_a - Y[28]) / R_ICAs;
    data->q_ICAr = (P_a - Y[29]) / R_ICAs;
    data->q_MCAl = (Y[28] - Y[22]) / (R_MCAs + (R_dml / 2.0));
    data->q_MCAr = (Y[29] - Y[25]) / (R_MCAs + (R_dmr / 2.0));
    data->q_ACAl = (P_ACA1l - Y[23]) / (R_ACA2s + (R_dal / 2.0));
    data->q_ACAr = (P_ACA1r - Y[26]) / (R_ACA2s + (R_dar / 2.0));
    data->q_PCAl = (P_PCA1l - Y[24]) / (R_PCA2s + (R_dpl / 2.0));
    data->q_PCAr = (P_PCA1r - Y[27]) / (R_PCA2s + (R_dpr / 2.0));

    // double data->q_PCoAl, data->q_PCoAr;
    if (PCoAl == 0) {
        data->q_PCoAl = (P_PCA1l - Y[28]) / R_PCoAs;
    } else if (PCoAl == 1) {
        data->q_PCoAl = 0.0;
    }

    if (PCoAr == 0) {
        data->q_PCoAr = (P_PCA1r - Y[29]) / R_PCoAs;
    } else if (PCoAr == 1) {
        data->q_PCoAr = 0.0;
    }

    // eq 5 V_dms
    double dV_dml = (Y[28] - Y[22]) / (R_MCAs + (R_dml / 2.0)) - (Y[22] - Y[20]) / (R_dml / 2.0) + (Y[23] - Y[22]) / R_cams + (Y[24] - Y[22]) / R_cpms; // for l
    double dV_dmr = (Y[29] - Y[25]) / (R_MCAs + (R_dmr / 2.0)) - (Y[25] - Y[20]) / (R_dmr / 2.0) + (Y[26] - Y[25]) / R_cams + (Y[27] - Y[25]) / R_cpms; // for r

    // eq 6 V_das
    double dV_dal = (P_ACA1l - Y[23]) / (R_ACA2s + (R_dal / 2.0)) - (Y[23] - Y[20]) / (R_dal / 2.0) + (Y[22] - Y[23]) / R_cams + (Y[26] - Y[23]) / R_caa; // for l
    double dV_dar = (P_ACA1r - Y[26]) / (R_ACA2s + (R_dar / 2.0)) - (Y[26] - Y[20]) / (R_dar / 2.0) + (Y[25] - Y[26]) / R_cams + (Y[23] - Y[26]) / R_caa; // for r

    // eq 7 V_dps
    double dV_dpl = (P_PCA1l - Y[24]) / (R_PCA2s + (R_dpl / 2.0)) - (Y[24] - Y[20]) / (R_dpl / 2.0) + (Y[22] - Y[24]) / R_cpms + (Y[27] - Y[24]) / R_cpp;
    double dV_dpr = (P_PCA1r - Y[27]) / (R_PCA2s + (R_dpr / 2.0)) - (Y[27] - Y[20]) / (R_dpr / 2.0) + (Y[25] - Y[27]) / R_cpms + (Y[24] - Y[27]) / R_cpp;

    // eq 4
    double C_ic = 1.0 / (k_E * Y[19]);

    // eq 1: dP_ic
    dY[19] = (dV_dml + dV_dpl + dV_dal + dV_dmr + dV_dpr + dV_dar + dV_vi + qf - qo) / C_ic;

    // double P_c = (Y[22]/(R_dml/2.0) + Y[23]/(R_dal/2.0) + Y[24]/(R_dpl/2.0) + Y[25]/(R_dmr/2.0) + Y[26]/(R_dar/2.0) + Y[27]/(R_dpr/2.0) + Y[21]/R_pv - qf)/(1.0/(R_dml/2.0) + 1.0/(R_dal/2.0) + 1.0/(R_dpl/2.0) + 1.0/(R_dmr/2.0) + 1.0/(R_dar/2.0) + 1.0/(R_dpr/2.0) + 1.0/R_pv);

    data->q_ml = (Y[22] - Y[20]) / (R_dml / 2.0);
    data->q_al = (Y[23] - Y[20]) / (R_dal / 2.0);
    data->q_pl = (Y[24] - Y[20]) / (R_dpl / 2.0);
    data->q_mr = (Y[25] - Y[20]) / (R_dmr / 2.0);
    data->q_ar = (Y[26] - Y[20]) / (R_dar / 2.0);
    data->q_pr = (Y[27] - Y[20]) / (R_dpr / 2.0);

    // dP_c
    double P_c = (data->q_ml + data->q_al + data->q_pl + data->q_mr + data->q_ar + data->q_pr - qf) * R_pv + Y[21];
    dY[20] = (P_c - Y[20]) / DELTAT;
    // dY[20] = (data->q_ml + data->q_al + data->q_pl + data->q_mr + data->q_ar + data->q_pr - qf - (Y[20] - Y[21])/R_pv)/0.1;

    // eq. 12
    double C_vi = 1.0 / (k_ven * (Y[21] - Y[19] - P_v1));

    // eq. 11 (reformatted): dP_v
    dY[21] = dV_vi / C_vi + dY[19];

    // eq 13 (reformatted): dP_djs
    dY[22] = dV_dml / Y[43] + dY[19] - dY[43] * (Y[22] - Y[19]) / Y[43];
    dY[23] = dV_dal / Y[44] + dY[19] - dY[44] * (Y[23] - Y[19]) / Y[44];
    dY[24] = dV_dpl / Y[45] + dY[19] - dY[45] * (Y[24] - Y[19]) / Y[45];
    dY[25] = dV_dmr / Y[46] + dY[19] - dY[46] * (Y[25] - Y[19]) / Y[46];
    dY[26] = dV_dar / Y[47] + dY[19] - dY[47] * (Y[26] - Y[19]) / Y[47];
    dY[27] = dV_dpr / Y[48] + dY[19] - dY[48] * (Y[27] - Y[19]) / Y[48];

    // eq 14 (reformatted): dP_ICAl, dP_ICAr
    if (A1l == 0) {
        dY[28] = (data->q_ICAl + data->q_PCoAl + (Y[22] - Y[28]) / (R_MCAs + R_dml / 2.0) + (P_ACA1l - Y[28]) / R_ACA1s) / C_ICAs + dY[19];
    } else if (A1l == 1) {
        dY[28] = (data->q_ICAl + data->q_PCoAl + (Y[22] - Y[28]) / (R_MCAs + R_dml / 2.0)) / C_ICAs + dY[19];
    }
    dY[29] = (data->q_ICAr + data->q_PCoAr + (Y[25] - Y[29]) / (R_MCAs + R_dmr / 2.0) + (P_ACA1r - Y[29]) / R_ACA1s) / C_ICAs + dY[19];

    // eq 15: dP_BA
    if (P1l == 0) {
        dY[30] = ((P_a - Y[30]) / R_BA + (P_PCA1l - Y[30]) / R_PCA1s + (P_PCA1r - Y[30]) / R_PCA1s) / C_BA;
    } else if (P1l == 1) {
        dY[30] = ((P_a - Y[30]) / R_BA + (P_PCA1r - Y[30]) / R_PCA1s) / C_BA;
    }

    // eq 18
    // double CBF;
    data->CBF = data->q_ml + data->q_al + data->q_pl + data->q_mr + data->q_ar + data->q_pr;

    // eq 19
    double V_MCAl = ((Y[28] - Y[22]) / (R_MCAs + R_dml / 2.0)) / (M_PI * pow(r_MCAl, 2));
    double V_MCAr = ((Y[29] - Y[25]) / (R_MCAs + R_dmr / 2.0)) / (M_PI * pow(r_MCAr, 2));

    // eq 21: dx_autjs
    dY[31] = (-Y[31] + G_autjs * ((data->q_ml - q_nms) / q_nms)) / tau_autjs;
    dY[32] = (-Y[32] + G_autjs * ((data->q_al - q_nas) / q_nas)) / tau_autjs;
    dY[33] = (-Y[33] + G_autjs * ((data->q_pl - q_nps) / q_nps)) / tau_autjs;
    dY[34] = (-Y[34] + G_autjs * ((data->q_mr - q_nms) / q_nms)) / tau_autjs;
    dY[35] = (-Y[35] + G_autjs * ((data->q_ar - q_nas) / q_nas)) / tau_autjs;
    dY[36] = (-Y[36] + G_autjs * ((data->q_pr - q_nps) / q_nps)) / tau_autjs;

    // eq 23
    double A_CO2ml = 1.0 / (1.0 + exp((-k_CO2 * (data->q_ml - q_nms) / q_nms) - b_CO2));
    double A_CO2al = 1.0 / (1.0 + exp((-k_CO2 * (data->q_al - q_nas) / q_nas) - b_CO2));
    double A_CO2pl = 1.0 / (1.0 + exp((-k_CO2 * (data->q_pl - q_nps) / q_nps) - b_CO2));
    double A_CO2mr = 1.0 / (1.0 + exp((-k_CO2 * (data->q_mr - q_nms) / q_nms) - b_CO2));
    double A_CO2ar = 1.0 / (1.0 + exp((-k_CO2 * (data->q_ar - q_nas) / q_nas) - b_CO2));
    double A_CO2pr = 1.0 / (1.0 + exp((-k_CO2 * (data->q_pr - q_nps) / q_nps) - b_CO2));

    // eq 22: dx_CO2js
    dY[37] = (-Y[37] + G_CO2js * A_CO2ml * log10(P_aCO2 / P_aCO2njs)) / tau_CO2js;
    dY[38] = (-Y[38] + G_CO2js * A_CO2al * log10(P_aCO2 / P_aCO2njs)) / tau_CO2js;
    dY[39] = (-Y[39] + G_CO2js * A_CO2pl * log10(P_aCO2 / P_aCO2njs)) / tau_CO2js;
    dY[40] = (-Y[40] + G_CO2js * A_CO2mr * log10(P_aCO2 / P_aCO2njs)) / tau_CO2js;
    dY[41] = (-Y[41] + G_CO2js * A_CO2ar * log10(P_aCO2 / P_aCO2njs)) / tau_CO2js;
    dY[42] = (-Y[42] + G_CO2js * A_CO2pr * log10(P_aCO2 / P_aCO2njs)) / tau_CO2js;

    // These are the flow in and out that connect the cerebral model to the body via the aorta and sup. vena cava
    double Q_Ceri = (P_a - Y[28]) / R_ICAs + (P_a - Y[29]) / R_ICAs + (P_a - Y[30]) / R_BA;
    double Q_Cero = (Y[21] - P_vs) / R_vs;
    data->temp = P_PCA1l;

    /***********************ODES***************************************************/
    /******************************************************************************/
    // %================== ALL ODEs ARE LISTED BELOW ==================

    /* Y[3] is Pup, pressure in upper body, head and neck.
                    eq. 22, appendix 2, p585.*/
    dY[0] = (data->Qupi - data->Qupo) / data->p_C[7];

    /*************** kidney ************************************/
    /*************** original formula **************************/
    dY[1] = (data->Qk1 - data->Qk2) / data->p_C[9]; // % Y[5-1] is Pk, Ck is non-zero. eq. 23, appendix 2, p585.

    /*************** detailed formula **************************/
    // for(i = 0; i<12; i++){
    //   dY[36+i] = (data->Q_kidneys[i] - data->Q_kidneys[12 + i])/(data->p_C[14+i]);
    // }
    //
    // // right and left feeding arteries to the kidneys.
    // double Q_kidRi_tot = 0.0;
    // double Q_kidLi_tot = 0.0;
    // for(i = 0; i<6; i++){
    //   Q_kidRi_tot = Q_kidRi_tot + data->Q_kidneys[i];
    //   Q_kidLi_tot = Q_kidLi_tot + data->Q_kidneys[6+i];
    // }
    // dY[48] = (data->QkRi - Q_kidRi_tot)/data->p_C[12]; // right.
    // dY[49] = (data->QkLi - Q_kidLi_tot)/data->p_C[13]; // left.
    /*************** END of  detailed kidney ************************************/

    // % Y[5] is Psp, splanchic pressure at inlet.
    dY[2] = (data->Qsp1 - data->Qsp2) / data->p_C[8];

    // % Y[6] is Pll, legs pressure.
    dY[3] = (data->Qll1 - data->Qll2) / data->p_C[10];

    // % Y[7] is Pab.
    dY[4] = (data->Qk2 + data->Qsp2 + data->Qll2 - data->Qab) / data->p_C[2];

    // % Y[8] is Pth: thoracic pressure. aka Y8 is Pbias. p1241, Heldt 2002 paper.
    dY[5] = (2.0 * M_PI * respRate) * cos(2.0 * M_PI * respRate * data->p_ursino[1]); // % This varies with respiratory variation, mean value taken from Heldt thesis 2.3 p.47 resprate units are (breath/s).

    /***************************** ODEs for Compliances **************************/
    dY[6] = (Clv - Y[6]) / DELTAT; // Left Ventricle Compliance
    dY[7] = (Crv - Y[7]) / DELTAT; // Right Ventricle Compliance
    dY[15] = (Cla - Y[15]) / DELTAT; // Left Atrial Compliance
    dY[16] = (Cra - Y[16]) / DELTAT; // Right Atrial Compliance

    /****************************************************************************************************************************/
    // % Y[11] is LV pressure.
    dY[8] = (1.0 / Y[6]) * ((Y[5] - Y[8]) * dY[6] + data->Qlao - data->Qlo) + dY[5];

    // Y34 is left atrial pressure.
    dY[17] = (1.0 / Y[15]) * ((Y[5] - Y[17]) * dY[15] + data->Qli - data->Qlao) + dY[5];

    // Y[15] is Right Ventricle Pressure
    dY[12] = (1.0 / Y[7]) * ((Y[5] - Y[12]) * dY[7] + data->Qrao - data->Qro) + dY[5];

    // Y[35] is right atrial pressure
    dY[18] = (1.0 / Y[16]) * ((Y[5] - Y[18]) * dY[16] + data->Qinf + data->Qsup - data->Qrao) + dY[5];

    // % Y[12] is Pa, ascending aorta pressure.
    // if(data->argv[3] == 1.0){
    //   dY[12]    = (1.0/data->p_C[6]) * (data->Qlo + ( -(Fa - Rv) - QF + Qinfused) - (data->Q_ao_bra + data->Q_ao_tho)) + 0.333*dY[8];
    // }else if(data->argv[3] == 0.0){
    //   dY[12]    = (1.0/data->p_C[6]) * (data->Qlo - (data->Q_ao_bra + data->Q_ao_tho)) + 0.333*dY[8];
    // }
    // Y[12] is systemic arterial pressure (aorta)
    dY[9] = (1.0 / data->p_C[6]) * (data->Qlo - (data->Qupi + data->Qsp1 + data->Qk1 + data->Qll1 + Q_Ceri)) + 0.333 * dY[5];

    // % Y[13] is Psup, eq. 27 appendix 2.
    dY[10] = ((data->Qupo + Q_Cero) - data->Qsup) / data->p_C[1] + dY[5];
    // % Y[14] is Pinf, eq. 28 appendix 2.
    dY[11] = (data->Qab - data->Qinf) / data->p_C[3] + dY[5];

    // % Y[16] is Ppa. Pulmonary circulation inlet pressure.
    dY[13] = (data->Qro - data->Qpa) / data->p_C[4] + dY[5]; // % Y[17-1] is Ppa, eq. 30 appendix 2.

    // % Y[17] is Ppv. Left ventricle inlet pressure.
    dY[14] = (data->Qpa - data->Qli) / data->p_C[5] + dY[5]; // % Y[18-1] is Ppv, eq. 31 appendix 2.

    for (int i = 0; i < NEQ; i++)
        Ith(ydot, i + 1) = dY[i];
    return (0);
}
