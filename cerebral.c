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

double P_a = Y[9];
double P_aCO2 = data->p_ursino[5];
double P_vs = Y[10];

///////////////////////// CoW variations //////////////////////////////////////
// The following are settings that determine which blood vessels are to be excluded from the model
int A1l, P1l, PCoAl, PCoAr;

if (data->CoW == 0){
  // Normal CoW
  A1l = 0;
  P1l = 0;
  PCoAl = 0;
  PCoAr = 0;
} else if (data->CoW == 1){
  // Left PCoA missing; PCoA
  A1l = 0;
  P1l = 0;
  PCoAl = 1;
  PCoAr = 0;
} else if (data->CoW == 2){
  // Both PCoA missing; PCoAs
  A1l = 0;
  P1l = 0;
  PCoAl = 1;
  PCoAr = 1;
} else if (data->CoW == 3){
  // ACAl missing; A1
  A1l = 1;
  P1l = 0;
  PCoAl = 0;
  PCoAr = 0;
} else if (data->CoW == 4){
  // PCAl missing; P1
  A1l = 0;
  P1l = 1;
  PCoAl = 0;
  PCoAr = 0;
} else if (data->CoW == 5){
  //PCoA and contralateral PCA; PCoA, PCA
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

for (int i = 0; i<6; i++) {

// eq. 25
  if (Y[37+i] - Y[31+i] < 0.0){
      deltaC_d[i] = 2.0 * sat1;
  } else {
      deltaC_d[i] = 2.0 * sat2;
  }

  k_Cd[i] = C_d0[i] * deltaC_d[i] / 4.0;


  C_dqs[i] = C_d0[i]*((1.0 - deltaC_d[i]/2.0) + (1.0 + deltaC_d[i]/2.0) * exp((Y[37+i] - Y[31+i])/k_Cd[i])) / (1.0 + exp((Y[37+i] - Y[31+i])/k_Cd[i]));
  if (C_dqs[i] < 0.0001){
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
    R_vs = (Y[21] - P_vs) * R_vs1 /(Y[21] - Y[19]);
} else {
    R_vs = R_vs1;
}



// eq 2
double qf = 0.0;
if (Y[20] > Y[19]){
    qf = (Y[20] - Y[19])/ R_f;
} else {
    qf = 0.0;
}

// eq 3
double qo = 0.0;
if (Y[19] > P_vs) {
    qo = (Y[19] - P_vs)/ R_o;
} else {
    qo = 0.0;
}

// eq 9
double dV_vi = (Y[20] - Y[21])/R_pv - (Y[21] - P_vs)/R_vs + qo;


double C_d0ms = C_dn * W_MCAs;
double C_d0as = C_dn * W_ACA2s;
double C_d0ps = C_dn * W_PCA2s;

// eq 26
double R_dml = (k_Rms * pow(C_d0ms,2))/(pow(Y[43]*(Y[22] - Y[19]),2));
double R_dal = (k_Ras * pow(C_d0as,2))/(pow(Y[44]*(Y[23] - Y[19]),2));
double R_dpl = (k_Rps * pow(C_d0ps,2))/(pow(Y[45]*(Y[24] - Y[19]),2));
double R_dmr = (k_Rms * pow(C_d0ms,2))/(pow(Y[46]*(Y[25] - Y[19]),2));
double R_dar = (k_Ras * pow(C_d0as,2))/(pow(Y[47]*(Y[26] - Y[19]),2));
double R_dpr = (k_Rps * pow(C_d0ps,2))/(pow(Y[48]*(Y[27] - Y[19]),2));


// eq 27
// These values are either from poiseills formula OR R_dn / Weight factor
double R_PCA1s = 0.7640;//(8.0 * mu * l_PCA1ns)/(M_PI * pow(r_PCA1ns,4));
double R_PCA2s = 3.6063;//(92.5 - 58.5) / (q_n * W_PCA2s); //(8.0 * mu * l_PCA2ns)/(M_PI * pow(r_PCA2ns,4)); // R_dn * W_PCA2s;
double R_PCoAs = 90.9786;//(8.0 * mu * l_PCoAns)/(M_PI * pow(r_PCoAns,4));
double R_ACA1s = 3.7912;//(8.0 * mu * l_ACA1ns)/(M_PI * pow(r_ACA1ns,4));
double R_ACA2s = 1.6227;//(92.5 - 58.5) / (q_n * W_ACA2s); // (8.0 * mu * l_ACA2ns)/(M_PI * pow(r_ACA2ns,4));
double R_ACoA  = 14.9228;//(8.0 * mu * l_ACoAn)/(M_PI *pow(r_ACoAn,4));
double R_ICAs  = 0.5689;//(100.0 - 92.5)/(q_n * W_ICAs);   //(8.0 * mu * l_ICAns)/(M_PI * pow(r_ICAns,4)); // R_dn * W_ICAs;
double R_BA    = 0.4501;//(100.0 - 92.5)/(q_n * W_BA);   //(8.0 * mu * l_BAn)/(M_PI * pow(r_BAn,4)); // R_dn * W_BA;


double P_ICAns = 92.5; // Units: mmHg, ref. p. 962 last paragraph

// eq 20
double r_MCAl = r_MCAns * ((1.0/k_MCAs)*log(Y[28]/P_ICAns) + 1.0);
double r_MCAr = r_MCAns * ((1.0/k_MCAs)*log(Y[29]/P_ICAns) + 1.0);

double R_MCAs = 1.4419;//(8.0 * mu * l_MCAns)/(M_PI * pow(r_MCAl,4));

// eq 17 (system solved using matlab solve function)
double P_ACA1l = 0.0;
double P_ACA1r = 0.0;
if (A1l == 0){
  P_ACA1l = (4.0*Y[28]*R_ACoA*pow(R_ACA2s,2) + 4.0*Y[28]*R_ACA1s*pow(R_ACA2s,2) + 4.0*Y[29]*R_ACA1s*pow(R_ACA2s,2) + 4.0*Y[23]*R_ACoA*pow(R_ACA1s,2) + 4.0*Y[23]*pow(R_ACA1s,2)*R_ACA2s + 4.0*Y[26]*pow(R_ACA1s,2)*R_ACA2s + 2.0*Y[23]*pow(R_ACA1s,2)*R_dar + 2.0*Y[26]*pow(R_ACA1s,2)*R_dal + 4.0*Y[28]*R_ACoA*R_ACA1s*R_ACA2s + 4.0*Y[23]*R_ACoA*R_ACA1s*R_ACA2s + 2.0*Y[28]*R_ACoA*R_ACA1s*R_dal + 2.0*Y[28]*R_ACoA*R_ACA2s*R_dal + 2.0*Y[28]*R_ACA1s*R_ACA2s*R_dal + 2.0*Y[28]*R_ACoA*R_ACA2s*R_dar + 2.0*Y[28]*R_ACA1s*R_ACA2s*R_dar + 2.0*Y[29]*R_ACA1s*R_ACA2s*R_dal + 2.0*Y[29]*R_ACA1s*R_ACA2s*R_dar + 2.0*Y[23]*R_ACoA*R_ACA1s*R_dar + Y[28]*R_ACoA*R_dal*R_dar + Y[28]*R_ACA1s*R_dal*R_dar + Y[29]*R_ACA1s*R_dal*R_dar)/(4.0*R_ACoA*pow(R_ACA1s,2) + 4.0*R_ACoA*pow(R_ACA2s,2) + 8.0*R_ACA1s*pow(R_ACA2s,2) + 8.0*pow(R_ACA1s,2)*R_ACA2s + 2.0*pow(R_ACA1s,2)*R_dal + 2.0*pow(R_ACA1s,2)*R_dar + 8.0*R_ACoA*R_ACA1s*R_ACA2s + 2.0*R_ACoA*R_ACA1s*R_dal + 2.0*R_ACoA*R_ACA2s*R_dal + 2.0*R_ACoA*R_ACA1s*R_dar + 4.0*R_ACA1s*R_ACA2s*R_dal + 2.0*R_ACoA*R_ACA2s*R_dar + 4.0*R_ACA1s*R_ACA2s*R_dar + R_ACoA*R_dal*R_dar + 2.0*R_ACA1s*R_dal*R_dar);

  P_ACA1r = (4.0*Y[28]*R_ACA1s*pow(R_ACA2s,2) + 4.0*Y[29]*R_ACoA*pow(R_ACA2s,2) + 4.0*Y[29]*R_ACA1s*pow(R_ACA2s,2) + 4.0*Y[23]*pow(R_ACA1s,2)*R_ACA2s + 4.0*Y[26]*R_ACoA*pow(R_ACA1s,2) + 4.0*Y[26]*pow(R_ACA1s,2)*R_ACA2s + 2.0*Y[23]*pow(R_ACA1s,2)*R_dar + 2.0*Y[26]*pow(R_ACA1s,2)*R_dal + 4.0*Y[29]*R_ACoA*R_ACA1s*R_ACA2s + 4.0*Y[26]*R_ACoA*R_ACA1s*R_ACA2s + 2.0*Y[28]*R_ACA1s*R_ACA2s*R_dal + 2.0*Y[29]*R_ACoA*R_ACA2s*R_dal + 2.0*Y[28]*R_ACA1s*R_ACA2s*R_dar + 2.0*Y[29]*R_ACoA*R_ACA1s*R_dar + 2.0*Y[29]*R_ACA1s*R_ACA2s*R_dal + 2.0*Y[29]*R_ACoA*R_ACA2s*R_dar + 2.0*Y[29]*R_ACA1s*R_ACA2s*R_dar + 2.0*Y[26]*R_ACoA*R_ACA1s*R_dal + Y[28]*R_ACA1s*R_dal*R_dar + Y[29]*R_ACoA*R_dal*R_dar + Y[29]*R_ACA1s*R_dal*R_dar)/(4.0*R_ACoA*pow(R_ACA1s,2) + 4.0*R_ACoA*pow(R_ACA2s,2) + 8.0*R_ACA1s*pow(R_ACA2s,2) + 8.0*pow(R_ACA1s,2)*R_ACA2s + 2.0*pow(R_ACA1s,2)*R_dal + 2.0*pow(R_ACA1s,2)*R_dar + 8.0*R_ACoA*R_ACA1s*R_ACA2s + 2.0*R_ACoA*R_ACA1s*R_dal + 2.0*R_ACoA*R_ACA2s*R_dal + 2.0*R_ACoA*R_ACA1s*R_dar + 4.0*R_ACA1s*R_ACA2s*R_dal + 2.0*R_ACoA*R_ACA2s*R_dar + 4.0*R_ACA1s*R_ACA2s*R_dar + R_ACoA*R_dal*R_dar + 2.0*R_ACA1s*R_dal*R_dar);

} else if (A1l == 1) {
  P_ACA1l = (4.0*Y[29]*pow(R_ACA2s,2) + 4.0*Y[23]*R_ACoA*R_ACA1s + 4.0*Y[23]*R_ACoA*R_ACA2s + 4.0*Y[23]*R_ACA1s*R_ACA2s + 4.0*Y[26]*R_ACA1s*R_ACA2s + 2.0*Y[29]*R_ACA2s*R_dal + 2.0*Y[29]*R_ACA2s*R_dar + 2.0*Y[23]*R_ACoA*R_dar + 2.0*Y[23]*R_ACA1s*R_dar + 2.0*Y[26]*R_ACA1s*R_dal + Y[29]*R_dal*R_dar)/(4.0*R_ACoA*R_ACA1s + 4.0*R_ACoA*R_ACA2s + 8.0*R_ACA1s*R_ACA2s + 2.0*R_ACA1s*R_dal + 2.0*R_ACoA*R_dar + 2.0*R_ACA2s*R_dal + 2.0*R_ACA1s*R_dar + 2.0*R_ACA2s*R_dar + R_dal*R_dar + 4.0*pow(R_ACA2s,2));

  P_ACA1r = (4.0*Y[29]*pow(R_ACA2s,2) + 4.0*Y[29]*R_ACoA*R_ACA2s + 4.0*Y[23]*R_ACA1s*R_ACA2s + 4.0*Y[26]*R_ACoA*R_ACA1s + 4.0*Y[26]*R_ACA1s*R_ACA2s + 2.0*Y[29]*R_ACoA*R_dar + 2.0*Y[29]*R_ACA2s*R_dal + 2.0*Y[29]*R_ACA2s*R_dar + 2.0*Y[23]*R_ACA1s*R_dar + 2.0*Y[26]*R_ACA1s*R_dal + Y[29]*R_dal*R_dar)/(4.0*R_ACoA*R_ACA1s + 4.0*R_ACoA*R_ACA2s + 8.0*R_ACA1s*R_ACA2s + 2.0*R_ACA1s*R_dal + 2.0*R_ACoA*R_dar + 2.0*R_ACA2s*R_dal + 2.0*R_ACA1s*R_dar + 2.0*R_ACA2s*R_dar + R_dal*R_dar + 4.0*pow(R_ACA2s,2));

}


// P_ACA1l = (P_dal/(R_ACA2s + R_dal/2) + P_ICAl/R_ACA1s + (P_dar/(R_ACA2s + R_dar/2) + P_ICAr/R_ACA1s) / (1/(R_ACA2s + R_dar/2) + 1/R_ACA1s + 1/R_ACoA)) / ((1/(R_ACA2s + R_dal/2) + 1/R_ACA1s + 1/R_ACoA) - (1/R_ACoA) * (1/(R_ACA2s + R_dar/2) + 1/R_ACA1s + 1/R_ACoA));
// P_ACA1r = (P_dar/(R_ACA2s + R_dar/2) + P_ICAr/R_ACA1s + (P_dal/(R_ACA2s + R_dal/2) + P_ICAl/R_ACA1s) / (1/(R_ACA2s + R_dal/2) + 1/R_ACA1s + 1/R_ACoA)) / ((1/(R_ACA2s + R_dar/2) + 1/R_ACA1s + 1/R_ACoA) - (1/R_ACoA) * (1/(R_ACA2s + R_dal/2) + 1/R_ACA1s + 1/R_ACoA));
//
// double P_ACA1l = (Y[23]/(R_ACA2s + R_dal/2.0) + Y[28]/R_ACA1s + (Y[26]/(R_ACA2s + R_dar/2.0) + Y[29]/R_ACA1s)/(1.0/(R_ACA2s + R_dar/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA)) / (1.0/(R_ACA2s + R_dal/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA - 1.0/(R_ACoA*(1.0/(R_ACA2s + R_dar/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA)));
// double P_ACA1r = (Y[26]/(R_ACA2s + R_dar/2.0) + Y[29]/R_ACA1s + (Y[23]/(R_ACA2s + R_dal/2.0) + Y[28]/R_ACA1s)/(1.0/(R_ACA2s + R_dal/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA)) / (1.0/(R_ACA2s + R_dar/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA - 1.0/(R_ACoA*(1.0/(R_ACA2s + R_dal/2.0) + 1.0/R_ACA1s + 1.0/R_ACoA)));


// eq 16
double P_PCA1l = 0.0;
if (P1l == 0){
  if (PCoAl == 0){
    P_PCA1l = (Y[28]/R_PCoAs + Y[30]/R_PCA1s + Y[24]/(R_PCA2s + R_dpl/2.0))/(1.0/R_PCoAs + 1.0/R_PCA1s + 1.0/(R_PCA2s + R_dpl/2.0));
  } else if (PCoAl == 1){
    P_PCA1l = (Y[30]/R_PCA1s + Y[24]/(R_PCA2s + R_dpl/2.0))/(1.0/R_PCA1s + 1.0/(R_PCA2s + R_dpl/2.0));
  }
} else if (P1l == 1) {
  if (PCoAl == 0){
    P_PCA1l = (Y[28]/R_PCoAs + Y[24]/(R_PCA2s + R_dpl/2.0))/(1.0/R_PCoAs + 1.0/(R_PCA2s + R_dpl/2.0));
  } else if (PCoAl == 1){
    P_PCA1l = (Y[24]/(R_PCA2s + R_dpl/2.0))/(1.0/(R_PCA2s + R_dpl/2.0));
  }
}
double P_PCA1r;
if (PCoAr == 0){
  P_PCA1r = (Y[29]/R_PCoAs + Y[30]/R_PCA1s + Y[27]/(R_PCA2s + R_dpr/2.0))/(1.0/R_PCoAs + 1.0/R_PCA1s + 1.0/(R_PCA2s + R_dpr/2.0));
} else if (PCoAr == 1){
  P_PCA1r = (Y[30]/R_PCA1s + Y[27]/(R_PCA2s + R_dpr/2.0))/(1.0/R_PCA1s + 1.0/(R_PCA2s + R_dpr/2.0));
}

// Flows
data->q_ICAl = (P_a - Y[28])/R_ICAs;
data->q_ICAr = (P_a - Y[29])/R_ICAs;
data->q_MCAl = (Y[28] - Y[22])/(R_MCAs + (R_dml/2.0));
data->q_MCAr = (Y[29] - Y[25])/(R_MCAs + (R_dmr/2.0));
data->q_ACAl = (P_ACA1l - Y[23])/(R_ACA2s + (R_dal/2.0));
data->q_ACAr = (P_ACA1r - Y[26])/(R_ACA2s + (R_dar/2.0));
data->q_PCAl = (P_PCA1l - Y[24])/(R_PCA2s + (R_dpl/2.0));
data->q_PCAr = (P_PCA1r - Y[27])/(R_PCA2s + (R_dpr/2.0));

// double data->q_PCoAl, data->q_PCoAr;
if (PCoAl == 0){
  data->q_PCoAl = (P_PCA1l - Y[28])/R_PCoAs;
} else if (PCoAl == 1){
  data->q_PCoAl = 0.0;
}

if (PCoAr == 0){
  data->q_PCoAr = (P_PCA1r - Y[29])/R_PCoAs;
} else if (PCoAr == 1){
  data->q_PCoAr = 0.0;
}

// eq 5 V_dms
double dV_dml = (Y[28] - Y[22])/(R_MCAs + (R_dml/2.0)) - (Y[22] - Y[20])/(R_dml/2.0) + (Y[23] - Y[22])/R_cams + (Y[24] - Y[22])/R_cpms; // for l
double dV_dmr = (Y[29] - Y[25])/(R_MCAs + (R_dmr/2.0)) - (Y[25] - Y[20])/(R_dmr/2.0) + (Y[26] - Y[25])/R_cams + (Y[27] - Y[25])/R_cpms; // for r

// eq 6 V_das
double dV_dal = (P_ACA1l - Y[23])/(R_ACA2s + (R_dal/2.0)) - (Y[23] - Y[20])/(R_dal/2.0) + (Y[22] - Y[23])/R_cams + (Y[26] - Y[23])/R_caa; // for l
double dV_dar = (P_ACA1r - Y[26])/(R_ACA2s + (R_dar/2.0)) - (Y[26] - Y[20])/(R_dar/2.0) + (Y[25] - Y[26])/R_cams + (Y[23] - Y[26])/R_caa; // for r

// eq 7 V_dps
double dV_dpl = (P_PCA1l - Y[24])/(R_PCA2s + (R_dpl/2.0)) - (Y[24] - Y[20])/(R_dpl/2.0) + (Y[22] - Y[24])/R_cpms + (Y[27] - Y[24])/R_cpp;
double dV_dpr = (P_PCA1r - Y[27])/(R_PCA2s + (R_dpr/2.0)) - (Y[27] - Y[20])/(R_dpr/2.0) + (Y[25] - Y[27])/R_cpms + (Y[24] - Y[27])/R_cpp;

// eq 4
double C_ic = 1.0/(k_E * Y[19]);

// eq 1: dP_ic
dY[19] = (dV_dml + dV_dpl + dV_dal + dV_dmr + dV_dpr + dV_dar + dV_vi + qf - qo)/C_ic;


// double P_c = (Y[22]/(R_dml/2.0) + Y[23]/(R_dal/2.0) + Y[24]/(R_dpl/2.0) + Y[25]/(R_dmr/2.0) + Y[26]/(R_dar/2.0) + Y[27]/(R_dpr/2.0) + Y[21]/R_pv - qf)/(1.0/(R_dml/2.0) + 1.0/(R_dal/2.0) + 1.0/(R_dpl/2.0) + 1.0/(R_dmr/2.0) + 1.0/(R_dar/2.0) + 1.0/(R_dpr/2.0) + 1.0/R_pv);

data->q_ml = (Y[22] - Y[20])/(R_dml/2.0);
data->q_al = (Y[23] - Y[20])/(R_dal/2.0);
data->q_pl = (Y[24] - Y[20])/(R_dpl/2.0);
data->q_mr = (Y[25] - Y[20])/(R_dmr/2.0);
data->q_ar = (Y[26] - Y[20])/(R_dar/2.0);
data->q_pr = (Y[27] - Y[20])/(R_dpr/2.0);


// dP_c
double P_c = (data->q_ml + data->q_al + data->q_pl + data->q_mr + data->q_ar + data->q_pr - qf) * R_pv + Y[21];
dY[20] = (P_c - Y[20])/DELTAT;
// dY[20] = (data->q_ml + data->q_al + data->q_pl + data->q_mr + data->q_ar + data->q_pr - qf - (Y[20] - Y[21])/R_pv)/0.1;

// eq. 12
double C_vi = 1.0/(k_ven * (Y[21] - Y[19] - P_v1));

// eq. 11 (reformatted): dP_v
dY[21] = dV_vi/C_vi + dY[19];


// eq 13 (reformatted): dP_djs
dY[22] = dV_dml/Y[43] + dY[19] - dY[43]*(Y[22] - Y[19])/Y[43];
dY[23] = dV_dal/Y[44] + dY[19] - dY[44]*(Y[23] - Y[19])/Y[44];
dY[24] = dV_dpl/Y[45] + dY[19] - dY[45]*(Y[24] - Y[19])/Y[45];
dY[25] = dV_dmr/Y[46] + dY[19] - dY[46]*(Y[25] - Y[19])/Y[46];
dY[26] = dV_dar/Y[47] + dY[19] - dY[47]*(Y[26] - Y[19])/Y[47];
dY[27] = dV_dpr/Y[48] + dY[19] - dY[48]*(Y[27] - Y[19])/Y[48];

// eq 14 (reformatted): dP_ICAl, dP_ICAr
if (A1l == 0){
  dY[28] = (data->q_ICAl + data->q_PCoAl + (Y[22] - Y[28])/(R_MCAs + R_dml/2.0) + (P_ACA1l - Y[28])/R_ACA1s)/C_ICAs + dY[19];
} else if (A1l == 1) {
  dY[28] = (data->q_ICAl + data->q_PCoAl + (Y[22] - Y[28])/(R_MCAs + R_dml/2.0))/C_ICAs + dY[19];
}
dY[29] = (data->q_ICAr + data->q_PCoAr + (Y[25] - Y[29])/(R_MCAs + R_dmr/2.0) + (P_ACA1r - Y[29])/R_ACA1s)/C_ICAs + dY[19];

// eq 15: dP_BA
if (P1l == 0){
  dY[30] = ((P_a - Y[30])/R_BA + (P_PCA1l - Y[30])/R_PCA1s + (P_PCA1r - Y[30])/R_PCA1s)/C_BA;
} else if (P1l == 1) {
  dY[30] = ((P_a - Y[30])/R_BA + (P_PCA1r - Y[30])/R_PCA1s)/C_BA;
}

// eq 18
// double CBF;
data->CBF = data->q_ml + data->q_al + data->q_pl + data->q_mr + data->q_ar + data->q_pr;

// eq 19
double V_MCAl = ((Y[28] - Y[22])/(R_MCAs + R_dml/2.0))/(M_PI * pow(r_MCAl,2));
double V_MCAr = ((Y[29] - Y[25])/(R_MCAs + R_dmr/2.0))/(M_PI * pow(r_MCAr,2));

// eq 21: dx_autjs
dY[31] = (-Y[31] + G_autjs * ((data->q_ml - q_nms)/q_nms))/tau_autjs;
dY[32] = (-Y[32] + G_autjs * ((data->q_al - q_nas)/q_nas))/tau_autjs;
dY[33] = (-Y[33] + G_autjs * ((data->q_pl - q_nps)/q_nps))/tau_autjs;
dY[34] = (-Y[34] + G_autjs * ((data->q_mr - q_nms)/q_nms))/tau_autjs;
dY[35] = (-Y[35] + G_autjs * ((data->q_ar - q_nas)/q_nas))/tau_autjs;
dY[36] = (-Y[36] + G_autjs * ((data->q_pr - q_nps)/q_nps))/tau_autjs;

// eq 23
double A_CO2ml = 1.0/(1.0 + exp((-k_CO2 * (data->q_ml - q_nms)/q_nms) - b_CO2));
double A_CO2al = 1.0/(1.0 + exp((-k_CO2 * (data->q_al - q_nas)/q_nas) - b_CO2));
double A_CO2pl = 1.0/(1.0 + exp((-k_CO2 * (data->q_pl - q_nps)/q_nps) - b_CO2));
double A_CO2mr = 1.0/(1.0 + exp((-k_CO2 * (data->q_mr - q_nms)/q_nms) - b_CO2));
double A_CO2ar = 1.0/(1.0 + exp((-k_CO2 * (data->q_ar - q_nas)/q_nas) - b_CO2));
double A_CO2pr = 1.0/(1.0 + exp((-k_CO2 * (data->q_pr - q_nps)/q_nps) - b_CO2));

// eq 22: dx_CO2js
dY[37] = (-Y[37] + G_CO2js * A_CO2ml * log10(P_aCO2/P_aCO2njs))/ tau_CO2js;
dY[38] = (-Y[38] + G_CO2js * A_CO2al * log10(P_aCO2/P_aCO2njs))/ tau_CO2js;
dY[39] = (-Y[39] + G_CO2js * A_CO2pl * log10(P_aCO2/P_aCO2njs))/ tau_CO2js;
dY[40] = (-Y[40] + G_CO2js * A_CO2mr * log10(P_aCO2/P_aCO2njs))/ tau_CO2js;
dY[41] = (-Y[41] + G_CO2js * A_CO2ar * log10(P_aCO2/P_aCO2njs))/ tau_CO2js;
dY[42] = (-Y[42] + G_CO2js * A_CO2pr * log10(P_aCO2/P_aCO2njs))/ tau_CO2js;


// These are the flow in and out that connect the cerebral model to the body via the aorta and sup. vena cava
double Q_Ceri = (P_a - Y[28])/R_ICAs + (P_a - Y[29])/R_ICAs + (P_a - Y[30])/R_BA;
double Q_Cero = (Y[21] - P_vs)/R_vs;
data->temp = P_PCA1l;

// syms P_ACA1l P_ACA1r
//
// syms P_dal P_dar P_ICAr P_ICAl R_ACA2s R_ACA1s R_dal R_ACoA  R_dar
//
// eqns = [(P_ACA1l - P_dal)/(R_ACA2s + R_dal/2) == (P_ICAl - P_ACA1l)/R_ACA1s + (P_ACA1r - P_ACA1l)/R_ACoA, ...
//     (P_ACA1r - P_dar)/(R_ACA2s + R_dar/2) == (P_ICAr - P_ACA1r)/R_ACA1s + (P_ACA1l - P_ACA1r)/R_ACoA];
//
// S = solve(eqns, [P_ACA1l P_ACA1r]);
