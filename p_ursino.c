// our my parameters. *************************************************************************************************************
data->p_ursino[0] = 0.0; // you never use p0, but it must have a value, lest you confuse the program.
data->p_ursino[1]  = 0.0;            // time, seconds. data->p_Ursino[1] = 0.0;

data->p_ursino[2] = randomPar(12/70.0, &randomParIndex, parameterFile, randomPars); // 12/70.0;        // baseline respiration rate; units: s/breath; ref. Heldt 2002

data->p_ursino[3] = 0.0;  // this is loc_t passed back from RHS to the main function thru' userdata. loc_t is from 0 to p[11].


/*******************************************************************************
HEART CHAMBER COMPLIANCES
*******************************************************************************/
// units: ml/mmHg;
// Heldt 2002, page 1243, paragraph 1.
// JJJ TH NOV 3 Took these values from heldt appendix

data->Edias_lv  = randomPar(0.13, &randomParIndex, parameterFile, randomPars);//*atof(argv[4]);
data->Esys_lv   = randomPar(2.5, &randomParIndex, parameterFile, randomPars);//*atof(argv[5]);

data->Edias_rv  = randomPar(0.07, &randomParIndex, parameterFile, randomPars);//*atof(argv[6]);
data->Esys_rv   = randomPar(1.3, &randomParIndex, parameterFile, randomPars);//*atof(argv[7]);

data->Edias_la  = randomPar(0.5, &randomParIndex, parameterFile, randomPars);//*atof(argv[8]);
data->Esys_la   = randomPar(0.61, &randomParIndex, parameterFile, randomPars);//*atof(argv[9]);

data->Edias_ra  = randomPar(0.3, &randomParIndex, parameterFile, randomPars);//*atof(argv[10]);
data->Esys_ra   = randomPar(0.74, &randomParIndex, parameterFile, randomPars);//*atof(argv[11]);

if(atoi(argv[2])==1){ // switch on for AF.
  data->Esys_la   = data->Edias_la; // 0.61*atof(argv[11]);
  data->Esys_ra   = data->Edias_ra; // 0.74*atof(argv[13]);
}


/******************************************************************************
RESISTANCE PARAMETERS
******************************************************************************/
// units: (mmHg*s/ml);
// Table 2 of Heldt 2002.

// lowering of temperature increases resistances. p99 multiples all resistances. at 36C, the increase is approx. 14%.
// rescaling of low resistances caused integration to fail around t = 280.
// a bunch of resistances are time dependent.

// data->p_R[] = ; //
data->p_R[0]  = randomPar(0.0, &randomParIndex, parameterFile, randomPars); // Empty
data->p_R[1]  = randomPar(0.06, &randomParIndex, parameterFile, randomPars);         // Rsup:  superior vena cava;
data->p_R[2]  = randomPar(0.01, &randomParIndex, parameterFile, randomPars);                                // Rab:   abdominal vena cava;
data->p_R[3]  = randomPar(0.015, &randomParIndex, parameterFile, randomPars);                               // Rinf:  inferior vena cava;
data->p_R[4]  = randomPar( 0.005, &randomParIndex, parameterFile, randomPars);                             // Rao: right heart's atrio-ventricular valve (tricuspid valve)
data->p_R[5]  = randomPar( 0.003, &randomParIndex, parameterFile, randomPars);                             // Rro: Resistance of right heart outlet;
data->p_R[6]  = randomPar(0.08, &randomParIndex, parameterFile, randomPars);        // Rp: Resistance of pulmonary arteries;
data->p_R[7]  = randomPar(0.01, &randomParIndex, parameterFile, randomPars);                               // Rpv: Resistance of pulmonary veins;
data->p_R[8]  = randomPar(0.01, &randomParIndex, parameterFile, randomPars);                               // Rmv: Left heart's mitral value
data->p_R[9]  = randomPar(0.006, &randomParIndex, parameterFile, randomPars);                             // Rlo: Resistance of left heart outlet;

data->p_R[10] = randomPar(8.1, &randomParIndex, parameterFile, randomPars);//*atof(argv[16]);  //  Rup1: Resistance of upper body(1), originally 3.9 (Heldt); changed to 8.1
data->p_R[11] = randomPar(0.5, &randomParIndex, parameterFile, randomPars);                //  Rup2: Resistance of upper body (2), originally 0.23 (Heldt);
data->p_R[12] = randomPar(3.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[14]);  //  Rsp1: Resistance of splanchic circulation (1);
data->p_R[13] = randomPar(0.18, &randomParIndex, parameterFile, randomPars);                //  Rsp2: Resistance of splanchic circulation (2);
data->p_R[14] = randomPar(4.1, &randomParIndex, parameterFile, randomPars);//*atof(argv[12]);                  //  Changed to R_kidi  //R_Inlet_right Kidneys
data->p_R[15] = randomPar(0.3, &randomParIndex, parameterFile, randomPars);//*atof(argv[13]);                  //  Changed to R_kido  //_Inlet_left Kidneys

data->p_R[16] = randomPar(3.6, &randomParIndex, parameterFile, randomPars);//*atof(argv[15]);                // Rll1: Resistance of legs (1);
data->p_R[17] = randomPar(0.3, &randomParIndex, parameterFile, randomPars);                               // Rll2: Resistance of legs (2);

/*******************************************************************************
COMPARTMENT COMPLIANCES
*******************************************************************************/
// units: mL/mmHg;
// Heldt 2002 AJP. p. 1242.

data->p_C[0] = randomPar(0.0, &randomParIndex, parameterFile, randomPars);    // not used
data->p_C[1] = randomPar(15.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[25]);   // Csup: Capacitance of superior vena cava;
data->p_C[2] = randomPar(25.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[21]);   // Cab: Capacitance of abdominal veins;
data->p_C[3] = randomPar(2.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[24]);    // Cinf: Capacitance of inferior vena cava;
data->p_C[4] = randomPar(4.3, &randomParIndex, parameterFile, randomPars);//*atof(argv[26]);    // Cpa: Capacitance of pulmonary arteries;
data->p_C[5] = randomPar(8.4, &randomParIndex, parameterFile, randomPars);//*atof(argv[27]);    // Cpv: Capacitance of pulmonary veins;
data->p_C[6] = randomPar(2.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[22]);    // Ca: Capacitance of systemic artery, i.e. aorta;
data->p_C[7] = randomPar(7.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[23]);   // Cup: Capacitance of upper body, originally 8.0;
data->p_C[8] = randomPar(55.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[19]);  // Csp: Splanchnic capacitance;
data->p_C[9] = randomPar(15.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[17]);  // Changed to C_kid   //C_Inlet_right Kidneys
data->p_C[10] =randomPar(19.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[20]);   // Cll: Legs venous capacitance;

/*******************************************************************************
BAROREFLEX PARAMETERS
*******************************************************************************/

data->p_ursino[4]   = 0.0;	//t_cardCycleInit: initiation time of the current cardiac cycle

data->p_ursino[5]   = randomPar(40.0, &randomParIndex, parameterFile, randomPars); // Pbco2; units: mmHg
data->p_ursino[6]   = randomPar(87.0, &randomParIndex, parameterFile, randomPars); // Pb02; units: mmHg

// Baroreflex INTEGRATOR Parameters
data->p_ursino[7]  = 0.001;   // tau_aff: units: seconds
data->p_ursino[8]  = 1;       // G_aff
data->p_ursino[9]  = 0.001;   // tau_c: Same as tau_aff,  units: s
data->p_ursino[10]  = 0.0205;  // S_p
data->p_ursino[11]  = 6;       // PNA_max, units: Hz
data->p_ursino[12]  = 0.6;     // PNA_min, units: Hz
data->p_ursino[13]  = -0.0138; // S_s
data->p_ursino[14]  = 4.0;       // SNA_max, units: Hz
data->p_ursino[15]  = 1.12;    // SNA_min, units: Hz
data->p_ursino[16]  = 13.8;    // k1
data->p_ursino[17]  = 0.182;   // k2
data->p_ursino[18]  = 828;     // k3
data->p_ursino[19]  = 1.0;       // k4
data->p_ursino[20]  = -18.118; // k5

// Baroreflex EFFECTOR Parameters
data->p_ursino[21]  = 90.0; 			// G_k_s0, units: beats/min/Hz
data->p_ursino[22]  = 0.28; 		// k_k_s0
data->p_ursino[23]  = 3.0; 			// T_s, units: s
data->p_ursino[24]  = 60.0; 			// G_v0, units: beats/min/Hz 45
data->p_ursino[25]  = 0.4; 		// k_v0
data->p_ursino[26]  = 0.5; 		// T_v, units: s
data->p_ursino[27]  = 1.5; 		// tau_sigma_lv, units: s
data->p_ursino[28]  = 2; 			// T_e_lv, units: s
data->p_ursino[29]  = 0.45*1.5; 		// G_eff_lv, units: mmHg/ml/Hz

data->p_ursino[30]  = 1.5; 		// tau_sigma_rv, units: s
data->p_ursino[31]  = 2; 			// T_e_rv, units: s
data->p_ursino[32]  = 0.282*0.935; 	// G_eff_rv, units: mmHg/ml/Hz
data->p_ursino[33]  = 10.0; 			// tau_sigma_V, units: s
data->p_ursino[34]  = 5; 			// T_e_V, units: s
data->p_ursino[35]  = -275.0*13.0;//*G_factor; 		// G_eff_V, units: ml/Hz
data->p_ursino[36]  = 1.5; 		// tau_sigma_R, units: s
data->p_ursino[37]  = 3; 			// T_e_R, units: s
data->p_ursino[38]  = 0.2*0.55; 		// G_eff_R, units: mmHg/(ml/s)/Hz // EDIT G_eff_R increased from 0.2 to 0.21
data->p_ursino[39]  = 25.0; 			// tau_s, simplified from equation 23, units: s
data->p_ursino[40]  = 0.8; 		// tau_v, simplified from eq. 28; units: s


data->p_ursino[41] 	= randomPar(75.0, &randomParIndex, parameterFile, randomPars);//*atof(argv[29]); 		// HR0 units: bpm

data->p_ursino[42]  	= 6.0; 			// lambda; for AF calculation
if (atoi(argv[2]) < 1){ // switch on for AF
  data->p_ursino[43]  	= 0.07; 			// coefficient of variation for HR
} else {
  data->p_ursino[43]  	= 0.24; 			// coefficient of variation for HR
}
// data->p_ursino[169]  	= 0.05; 			// coefficient of variation for HR
data->p_ursino[44]  	= 0.5; 			// tau; for AF calculation
data->p_ursino[45]  	= 0.1; 			// nu; for AF calculation

// This is a buffer for the past 2 RR intervals
data->RR[0]           = 0.8;
data->RR[1]           = 0.8;
data->RR[2]           = 0.8;

data->q_ml = 0.0;
data->q_al = 0.0;
data->q_pl = 0.0;
data->q_mr = 0.0;
data->q_ar = 0.0;
data->q_pr = 0.0;

data->q_ICAl = 0.0;
data->q_ICAr = 0.0;
data->q_MCAl = 0.0;
data->q_MCAr = 0.0;
data->q_ACAl = 0.0;
data->q_ACAr = 0.0;
data->q_PCAl = 0.0;
data->q_PCAr = 0.0;
data->q_PCoAl = 0.0;
data->q_PCoAr = 0.0;

data->temp = 0.0;


for(i = 0; i<(int)(5.0/DELTAT); i++){
  data->SNA_buffer[i] = 3.0;
}
for(i = 0; i<(int)(0.5/DELTAT); i++){
  data->PNA_buffer[i] = 3.0;
}


data->CoW = atoi(argv[3]);


// Parameters from table 1 of ref.
// Hemodynamic and hydrodynamic
data->P_cerebral[0] = randomPar(2.38 * pow(10,3), &randomParIndex, parameterFile, randomPars); // R_f, units mmHg s mL^-1
data->P_cerebral[1] = randomPar(526.3, &randomParIndex, parameterFile, randomPars); // R_o, units mmHg s mL^-1
data->P_cerebral[2] = randomPar(0.880, &randomParIndex, parameterFile, randomPars); // R_pv, units mmHg s mL^-1
data->P_cerebral[3] = randomPar(0.366, &randomParIndex, parameterFile, randomPars); // R_vs1, units mmHg s mL^-1
data->P_cerebral[4] = randomPar(120.0, &randomParIndex, parameterFile, randomPars); // R_cpms, units mmHg s mL^-1
data->P_cerebral[5] = randomPar(105.0, &randomParIndex, parameterFile, randomPars); // R_cams, units mmHg s mL^-1
data->P_cerebral[6] = randomPar(75.0, &randomParIndex, parameterFile, randomPars); // R_cpp,  mmHg s mL^-1
data->P_cerebral[7] = randomPar(22.0, &randomParIndex, parameterFile, randomPars); // R_caa, units mmHg s mL^-1
data->P_cerebral[8] = randomPar(3.4 * pow(10,-3), &randomParIndex, parameterFile, randomPars); // C_ICAs, units mL mmHg^-1
data->P_cerebral[9] = randomPar(1.7 * pow(10,-3), &randomParIndex, parameterFile, randomPars); // C_BA, units mL mmHg^-1
data->P_cerebral[10] = randomPar(0.205, &randomParIndex, parameterFile, randomPars); // r_ICAns, units cm
data->P_cerebral[11] = randomPar(0.17, &randomParIndex, parameterFile, randomPars); // r_BAn, units cm
data->P_cerebral[12] = randomPar(0.14, &randomParIndex, parameterFile, randomPars); // r_MCAns, units cm
data->P_cerebral[13] = randomPar(0.1, &randomParIndex, parameterFile, randomPars); // r_PCA1ns, units cm
data->P_cerebral[14] = randomPar(0.075, &randomParIndex, parameterFile, randomPars); // r_ACA1ns, units cm
data->P_cerebral[15] = randomPar(0.1, &randomParIndex, parameterFile, randomPars); // r_PCA2ns, units cm
data->P_cerebral[16] = randomPar(0.075, &randomParIndex, parameterFile, randomPars); // r_ACA2ns, units cm
data->P_cerebral[17] = randomPar(0.036, &randomParIndex, parameterFile, randomPars); // r_PCoAns, units cm
data->P_cerebral[18] = randomPar(0.04, &randomParIndex, parameterFile, randomPars); // r_ACoAn, units cm
data->P_cerebral[19] = randomPar(13.15, &randomParIndex, parameterFile, randomPars); // l_ICAns, units cm
data->P_cerebral[20] = randomPar(4.92, &randomParIndex, parameterFile, randomPars); // l_BAn, units cm
data->P_cerebral[21] = randomPar(7.25, &randomParIndex, parameterFile, randomPars); // l_MCAns, units cm
data->P_cerebral[22] = randomPar(1.0, &randomParIndex, parameterFile, randomPars); // l_PCA1ns, units cm
data->P_cerebral[23] = randomPar(1.57, &randomParIndex, parameterFile, randomPars); // l_ACA1ns, units cm
data->P_cerebral[24] = randomPar(4.72, &randomParIndex, parameterFile, randomPars); // l_PCA2ns, units cm
data->P_cerebral[25] = randomPar(0.672, &randomParIndex, parameterFile, randomPars); // l_ACA2ns, units cm
data->P_cerebral[26] = randomPar(2.0, &randomParIndex, parameterFile, randomPars); // l_PCoAns, units cm
data->P_cerebral[27] = randomPar(0.5, &randomParIndex, parameterFile, randomPars); // l_ACoAn, units cm
data->P_cerebral[28] = randomPar(0.155, &randomParIndex, parameterFile, randomPars); // k_ven, units mL^-1
data->P_cerebral[29] = randomPar(12.0, &randomParIndex, parameterFile, randomPars); // k_MCAs
data->P_cerebral[30] = randomPar(0.077, &randomParIndex, parameterFile, randomPars); // k_E, units mL^-1
data->P_cerebral[31] = randomPar(-2.5, &randomParIndex, parameterFile, randomPars); // P_v1, units mmHg
data->P_cerebral[32] = randomPar(10.1, &randomParIndex, parameterFile, randomPars); // V_dn, units mL
data->P_cerebral[33] = randomPar(200.0 * pow(10,-3), &randomParIndex, parameterFile, randomPars); // C_dn, units mL mmHg^-1
data->P_cerebral[34] = randomPar(5.4, &randomParIndex, parameterFile, randomPars); // R_dn, units mmHg s mL^-1
data->P_cerebral[35] = randomPar(13.1 * pow(10,3), &randomParIndex, parameterFile, randomPars); // k_R, units mmHg^-3 s mL^-1
data->P_cerebral[36] = randomPar(12.5, &randomParIndex, parameterFile, randomPars); // q_n, units mL s^-1
data->P_cerebral[37] = randomPar(0.38, &randomParIndex, parameterFile, randomPars); // W_ICAs
data->P_cerebral[38] = randomPar(0.24, &randomParIndex, parameterFile, randomPars); // W_BA
data->P_cerebral[39] = randomPar(0.12, &randomParIndex, parameterFile, randomPars); // W_PCA1s
data->P_cerebral[40] = randomPar(0.08, &randomParIndex, parameterFile, randomPars); // W_ACA1s
data->P_cerebral[41] = randomPar(0.30, &randomParIndex, parameterFile, randomPars); // W_MCAs
data->P_cerebral[42] = randomPar(0.12, &randomParIndex, parameterFile, randomPars); // W_PCA2s
data->P_cerebral[43] = randomPar(0.08, &randomParIndex, parameterFile, randomPars); // W_ACA2s

data->P_cerebral[44] = randomPar(0.004/133.322, &randomParIndex, parameterFile, randomPars); // mu, units poise

// Cerebrlovascular control mechanisms
data->P_cerebral[45] = randomPar(20.0, &randomParIndex, parameterFile, randomPars); // tau_autjs, units s
data->P_cerebral[46] = randomPar(0.9, &randomParIndex, parameterFile, randomPars); // G_autjs
data->P_cerebral[47] = randomPar(40.0, &randomParIndex, parameterFile, randomPars); // tau_CO2js, units s
data->P_cerebral[48] = randomPar(4.0, &randomParIndex, parameterFile, randomPars); // G_CO2js
data->P_cerebral[49] = randomPar(7.0, &randomParIndex, parameterFile, randomPars); // sat2
data->P_cerebral[50] = randomPar(0.4, &randomParIndex, parameterFile, randomPars); // sat1
data->P_cerebral[51] = randomPar(15.0, &randomParIndex, parameterFile, randomPars); // k_CO2
data->P_cerebral[52] = randomPar(0.5, &randomParIndex, parameterFile, randomPars); // b_CO2
data->P_cerebral[53] = randomPar(40.0, &randomParIndex, parameterFile, randomPars); // P_aCO2njs, units mmHg







for(i = 1; i<argc; i++){
  data->argv[i] = atof(argv[i]);
}
