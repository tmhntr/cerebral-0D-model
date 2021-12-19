// our my parameters. *************************************************************************************************************
data->p_ursino[0] = 0.0; // you never use p0, but it must have a value, lest you confuse the program.
data->p_ursino[1]  = 0.0;            // time, seconds. data->p_Ursino[1] = 0.0;

// data->p_ursino[11] = 60.0/70.0;      // Heart Rate; units (beats/s), baseline duration of cardiac cycle. Severi cycle length: 0.825 s per beat.
data->p_ursino[12] = 12/70.0;        // baseline respiration rate; units: s/breath; ref. Heldt 2002

/******************************************************************************/
//
// data->p_ursino[100] = 3 * pow(10,-5); // blood dynamic viscosity. Letcher et al. Am J Med. 70: 1195-1202. 1981. mmHg-s.
//
// // the 0.1 resistance is a placeholder - find justification when writing paper/doing simulations.
// data->p_ursino[101] = 0.1; // length of all microvasculature in cm, microvasculature is with resistance more than 0.1.
// data->p_ursino[102] = 10.0; // length of all large vessels in cm, microvasculature is with resistance less than 0.1.
//
// // the ascending aorta does not have a resistance as yet. For now, assign some constant value of radius to the aorta.
// // Literature value of ascending aorta diameter: 3 cm (Journal of the American College of Cardiology,  Volume 69, Issue 11 Supplement, March 2017, DOI: 10.1016/S0735-1097(17)35464-5; Paruchuri et al. Cardiology 2015;131:265-272).
// data->p_ursino[103] = 0.90; // ascending aorta diameter in cm, this is a large value as compared to what using lengths to work our radius gives in this model. The observed value is closer to 1.5 cm.
// data->p_ursino[104] = 0.007; // units PRU. Ascending aorta resistance according to table on p35 of Heldt thesis. Note that written in this way, it presumes certain things and will need to be generalised. This resistance is NOT feeding into the aorta flow yet. There is no aorta ODE yet. Written in this way, ascending aorta shear will depend on LV elastance.

data->p_ursino[105] = 0.0;  // this is loc_t passed back from RHS to the main function thru' userdata. loc_t is from 0 to p[11].


/*******************************************************************************
HEART CHAMBER COMPLIANCES
*******************************************************************************/
// units: ml/mmHg;
// Heldt 2002, page 1243, paragraph 1.
// JJJ TH NOV 3 Took these values from heldt appendix

data->Edias_lv  = 0.13;//*atof(argv[4]);
data->Esys_lv   = 2.5;//*atof(argv[5]);

data->Edias_rv  = 0.07;//*atof(argv[6]);
data->Esys_rv   = 1.3;//*atof(argv[7]);

data->Edias_la  = 0.5;//*atof(argv[8]);
data->Esys_la   = 0.61;//*atof(argv[9]);

data->Edias_ra  = 0.3;//*atof(argv[10]);
data->Esys_ra   = 0.74;//*atof(argv[11]);

if(atoi(argv[2])==1){ // switch on for AF.
  data->Edias_la  = 0.5;//*atof(argv[8]);
  data->Esys_la   = data->Edias_la; // 0.61*atof(argv[11]);

  data->Edias_ra  = 0.3;//*atof(argv[10]);
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
data->p_R[0]  = 0.0; // Empty
data->p_R[1]  = 0.06;         // Rsup:  superior vena cava;
data->p_R[2]  = 0.01;                                // Rab:   abdominal vena cava;
data->p_R[3]  = 0.015;                               // Rinf:  inferior vena cava;
data->p_R[4]  =  0.005;                             // Rao: right heart's atrio-ventricular valve (tricuspid valve)
data->p_R[5]  =  0.003;                             // Rro: Resistance of right heart outlet;
data->p_R[6]  = 0.08;        // Rp: Resistance of pulmonary arteries;
data->p_R[7]  = 0.01;                               // Rpv: Resistance of pulmonary veins;
data->p_R[8]  = 0.01;                               // Rmv: Left heart's mitral value
data->p_R[9]  = 0.006;                             // Rlo: Resistance of left heart outlet;
data->p_R[10] = 0.003;        // R_Brachiocephalic aorta
data->p_R[11] = 0.011;        // R_Thoracic Aorta
data->p_R[12] = 0.010;        // R_abdominal aorta

data->p_R[13] = 8.1;//*atof(argv[16]);  //  Rup1: Resistance of upper body(1), originally 3.9 (Heldt); changed to 8.1
data->p_R[14] = 0.5;                //  Rup2: Resistance of upper body (2), originally 0.23 (Heldt);
data->p_R[15] = 3.0;//*atof(argv[14]);  //  Rsp1: Resistance of splanchic circulation (1);
data->p_R[16] = 0.18;                //  Rsp2: Resistance of splanchic circulation (2);
data->p_R[17] = 4.1;//*atof(argv[12]);                  //  Changed to R_kidi  //R_Inlet_right Kidneys
data->p_R[18] = 0.3;//*atof(argv[13]);                  //  Changed to R_kido  //_Inlet_left Kidneys

for(i = 0; i<12; i++){
  data->p_R[19+i] = 4.1 * 12;   // R_kidi,1 to R_kidi,12
  }

for(i = 0; i<12; i++) data->p_R[31+i] = 0.3 * 12;  // R_kido,1 to R_kido,12

data->p_R[43] = 3.6;//*atof(argv[15]);                // Rll1: Resistance of legs (1);
data->p_R[44] = 0.3;                               // Rll2: Resistance of legs (2);

/*******************************************************************************
COMPARTMENT COMPLIANCES
*******************************************************************************/
// units: mL/mmHg;
// Heldt 2002 AJP. p. 1242.

data->p_C[0] = 0.0;    // not used
data->p_C[1] = 15.0;//*atof(argv[25]);   // Csup: Capacitance of superior vena cava;
data->p_C[2] = 25.0;//*atof(argv[21]);   // Cab: Capacitance of abdominal veins;
data->p_C[3] = 2.0;//*atof(argv[24]);    // Cinf: Capacitance of inferior vena cava;
data->p_C[4] = 4.3;//*atof(argv[26]);    // Cpa: Capacitance of pulmonary arteries;
data->p_C[5] = 8.4;//*atof(argv[27]);    // Cpv: Capacitance of pulmonary veins;
data->p_C[6] = 2.0;//*atof(argv[22]);    // Ca: Capacitance of systemic artery, i.e. aorta;
data->p_C[7] = 0.013;                  // C_Brachiocephalic aorta
data->p_C[8] = 0.021;                  // C_Thoracic Aorta
data->p_C[9] = 0.010;                  // C_abdominal aorta
data->p_C[10] = 7.0;//*atof(argv[23]);   // Cup: Capacitance of upper body, originally 8.0;
data->p_C[11] = 55.0;//*atof(argv[19]);  // Csp: Splanchnic capacitance;
data->p_C[12] = 15.0;//*atof(argv[17]);  // Changed to C_kid   //C_Inlet_right Kidneys
data->p_C[13] = 5.0;//*atof(argv[18]);  // C_Inlet_left Kidneys

for(i = 0; i<12; i++){                // Ckid1-12: Kidney capacitance;
  data->p_C[14+i] = (15.0/12.0);
}
data->p_C[26] = 19.0;//*atof(argv[20]);   // Cll: Legs venous capacitance;

/*******************************************************************************
BAROREFLEX PARAMETERS
*******************************************************************************/

data->p_ursino[106]   = 0.0;	//t_cardCycleInit: initiation time of the current cardiac cycle

data->p_ursino[107]   = 2.739; // SNA
data->p_ursino[108]   = 3.27; // PNA
data->p_ursino[109]   = 40.0; // Pbco2; units: mmHg
data->p_ursino[110]   = 87.0; // Pb02; units: mmHg
data->p_ursino[111]   = 0; // deltaHR_S
data->p_ursino[112]   = 0; // deltaHR_V
data->p_ursino[113]   = 1.15; // sigma_lv; units: mmHg/ml
data->p_ursino[114]   = 0.72; // sigma_rv; units: mmHg/ml
data->p_ursino[115]   = 9.3; // delta_sigma_V; units: ml
data->p_ursino[116]   = 0.7; // sigma_R; units: mmHg/(ml/s)

// ZPFV
data->p_ursino[117]   = 645;//650; // ZPFV_up: zero-pressure filling volume upper body veins; units: ml;
data->p_ursino[118]   = 30;//150; // ZPFV_kid;
data->p_ursino[119]   = 1146;//1300; // ZPFV_sp;
data->p_ursino[120]   = 716;//350; // ZPFV_ll
data->p_ursino[121]   = 79;//250; // ZPFV_ab
data->p_ursino[122]   = 33;//75; // ZPFV_inf
data->p_ursino[123]   = 16;//10; // ZPFV_sup

// Baroreflex INTEGRATOR Parameters
data->p_ursino[124]  = 0.001;   // tau_aff: units: seconds
data->p_ursino[125]  = 1;       // G_aff
data->p_ursino[126]  = 0.001;   // tau_c: Same as tau_aff,  units: s
data->p_ursino[127]  = 0.0205;  // S_p
data->p_ursino[128]  = 6;       // PNA_max, units: Hz
data->p_ursino[129]  = 0.6;     // PNA_min, units: Hz
data->p_ursino[130]  = -0.0138; // S_s
data->p_ursino[131]  = 4;       // SNA_max, units: Hz
data->p_ursino[132]  = 1.12;    // SNA_min, units: Hz
data->p_ursino[133]  = 13.8;    // k1
data->p_ursino[134]  = 0.182;   // k2
data->p_ursino[135]  = 828;     // k3
data->p_ursino[136]  = 1;       // k4
data->p_ursino[137]  = -18.118; // k5

// Baroreflex EFFECTOR Parameters
double G_factor = 1.0;//*atof(argv[28]);
data->p_ursino[138]  = 90.0; 			// G_k_s0, units: beats/min/Hz
data->p_ursino[139]  = 0.28; 		// k_k_s0
data->p_ursino[140]  = 3.0; 			// T_s, units: s
data->p_ursino[141]  = 60.0; 			// G_v0, units: beats/min/Hz 45
data->p_ursino[142]  = 0.4; 		// k_v0
data->p_ursino[143]  = 0.5; 		// T_v, units: s
data->p_ursino[144]  = 1.5; 		// tau_sigma_lv, units: s
data->p_ursino[145]  = 2; 			// T_e_lv, units: s
// All parameters with G_factor are baroreflex gains.
data->p_ursino[146]  = 0.45*1.5; 		// G_eff_lv, units: mmHg/ml/Hz

data->p_ursino[147]  = 1.5; 		// tau_sigma_rv, units: s
data->p_ursino[148]  = 2; 			// T_e_rv, units: s
data->p_ursino[149]  = 0.282*0.935; 	// G_eff_rv, units: mmHg/ml/Hz
data->p_ursino[150]  = 10.0; 			// tau_sigma_V, units: s
data->p_ursino[151]  = 5; 			// T_e_V, units: s
data->p_ursino[152]  = -275.0*13.0;//*G_factor; 		// G_eff_V, units: ml/Hz
data->p_ursino[153]  = 1.5; 		// tau_sigma_R, units: s
data->p_ursino[154]  = 3; 			// T_e_R, units: s
data->p_ursino[155]  = 0.2*0.55; 		// G_eff_R, units: mmHg/(ml/s)/Hz // EDIT G_eff_R increased from 0.2 to 0.21
data->p_ursino[156]  = 25.0; 			// tau_s, simplified from equation 23, units: s
data->p_ursino[157]  = 0.8; 		// tau_v, simplified from eq. 28; units: s
// state vectors for baroreflex
data->p_ursino[158]  	= 0.116; 			// x0_P_aff
data->p_ursino[159]  	= 0.0903; 			// x0_temp1
data->p_ursino[160] 	= 63.49; 			// x1_deltaHR_s;
data->p_ursino[161] 	= 2.92; 			// x1_deltaHR_v;
data->p_ursino[162] 	= 3.85; 			// x1_sigma_lv;
data->p_ursino[163] 	= 3.85; 			// x1_sigma_rv;
data->p_ursino[164] 	= 25.42; 			// x1_sigma_V;
data->p_ursino[165] 	= 3.84; 			// x1_sigma_R;
data->p_ursino[166] 	= -699.0; 			// sigma_V

data->p_ursino[167] 	= atof(argv[4]);//*atof(argv[29]); 		// HR0 units: bpm

data->p_ursino[168]  	= 6.0; 			// lambda; for AF calculation
if (atoi(argv[2]) < 1){ // switch on for AF
  data->p_ursino[169]  	= 0.07; 			// coefficient of variation for HR
} else {
  data->p_ursino[169]  	= 0.24; 			// coefficient of variation for HR
}
// data->p_ursino[169]  	= 0.05; 			// coefficient of variation for HR
data->p_ursino[170]  	= 0.5; 			// tau; for AF calculation
data->p_ursino[171]  	= 0.1; 			// nu; for AF calculation

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

for(i = 1; i<argc; i++){
  data->argv[i] = atof(argv[i]);
}
