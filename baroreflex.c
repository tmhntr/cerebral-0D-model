//  12 October, 2020
//  Implemented by Timothy Hunter
//  Lead Developer: Sanjay R. Kharche
//  Part of PM3
//  This code includes implementation of the baroreflex model from Lin et al.
//  2012.
//
//  Ref. DOI: 10.1177/0954411912451823
//  Constants are from Lin et al. 2012

// This model contains 8 ODEs


double MAP =  Y[9];                        // this is the aortic pressure, carotid pressures may be taken instead.
double Pbco2 =  data->p_ursino[90+19];      // this is the blood CO2 concentration
double Pbo2 =  data->p_ursino[91+19];       // this is the blood O2 concentration


//  Table 3 - Baroreflex control parameters
//  Afferent compartment:
double tau_aff =  data->p_ursino[105+19];   //  units: seconds
double G_aff   =  data->p_ursino[106+19];

//  Central Compartment:
double tau_c =  data->p_ursino[107+19];     //  Same as tau_aff,  units: s

//  Efferent compartment:
double S_p =  data->p_ursino[108+19];
double PNA_max =  data->p_ursino[109+19];   //  units: Hz
double PNA_min =  data->p_ursino[110+19];   //  units: Hz
double S_s =  data->p_ursino[111+19];
double SNA_max =  data->p_ursino[112+19];   //  units: Hz
double SNA_min =  data->p_ursino[113+19];   //  units: Hz

//  Table A3
double k1 =  data->p_ursino[114+19];
double k2 =  data->p_ursino[115+19];
double k3 =  data->p_ursino[116+19];
double k4 =  data->p_ursino[117+19];
double k5 =  data->p_ursino[118+19];

// P_aff
dY[49] = (-Y[49] + G_aff * (MAP))/tau_aff;


//  Central Compartment
double deltaMAP;
if ((Pbco2 > 40.0) && (Pbo2 < 104.0))
    deltaMAP = k1 + k2 * Pbco2 + k3 / Pbo2;
else if ((Pbco2 <= 40.0) && (Pbo2 < 104.0))
    deltaMAP = k1 + k2 * 40 + k3 / Pbo2;
else if ((Pbco2 > 40.0) && (Pbo2 >= 104.0))
    deltaMAP = k1 + k2 * Pbco2 + k3 / 104.0;


double P_demand = 90.0 + deltaMAP/100.0; //  Desrcibed in fig. 2 and on page 793

// for P_error
dY[50] = (-Y[50] + 1.0 * (P_demand))/tau_c; // second term in the P_error eq.

double P_error = - Y[49] + Y[50] ;//  eq. 13


//  Efferent Compartment:
//  eq. 31 chemoreceptor operating point
double deltaG_SNA;
if (Pbco2 > 40.0)
    deltaG_SNA = k4 * Pbco2 + k5;
else
    deltaG_SNA = k4 * 40.0 + k5;


double k_s = (SNA_max - SNA_min) / (4.0 * S_s); //  eq. 15
data->SNA = ((SNA_max + SNA_min * exp(P_error/k_s))/(1.0+ exp(P_error/k_s))) * (1.0 + deltaG_SNA/100.0); //  eq. 14

double S_v = S_p; //  paper uses inconsistent notation
double k_v = (PNA_max - PNA_min) / (4.0 * S_v); //  eq. 17
data->PNA = (PNA_max - PNA_min * exp(P_error/k_v)) / (1.0 + exp(P_error/k_v)); //  eq. 16


// Table 3 baroreflex control parameters
double G_k_s0 = data->p_ursino[119+19]; // units: beats/min/Hz
double k_k_s0 = data->p_ursino[120+19];
double T_s = data->p_ursino[121+19]; // units: s
double G_tau_s0 = 10.0; // beats/min/Hz
double k_tau_s0 = 0.2;


double G_v0 = data->p_ursino[122+19]; // units: beats/min/Hz 45
double k_v0 = data->p_ursino[123+19];
double T_v = data->p_ursino[124+19]; // units: s

// Myocardial constriction
double tau_sigma_lv = data->p_ursino[125+19]; // units: s
double T_e_lv = data->p_ursino[126+19]; // units: s
double G_eff_lv = data->p_ursino[127+19]; // units: mmHg/ml/Hz

double tau_sigma_rv = data->p_ursino[128+19]; // units: s
double T_e_rv = data->p_ursino[129+19]; // units: s
double G_eff_rv = data->p_ursino[130+19]; // units: mmHg/ml/Hz

// Veins
double tau_sigma_V = data->p_ursino[131+19]; // units: s
double T_e_V = data->p_ursino[132+19]; // units: s
double G_eff_V = data->p_ursino[133+19]; // units: ml/Hz

// Arterioles
double tau_sigma_R = data->p_ursino[134+19]; // units: s
double T_e_R = data->p_ursino[135+19]; // units: s
double G_eff_R = data->p_ursino[136+19]; // units: mmHg/(ml/s)/Hz // EDIT G_eff_R increased from 0.2 to 0.21

double tau_s;
if (data->SNA > data->SNA_buffer[0]){
  double k_tau_s = G_tau_s0 / k_tau_s0;
  tau_s = 10.0; //G_tau_s0 * exp(-k_tau_s * data->SNA); this is a substituted value based on G_tau_s0 as given eqation goes to zero and crashes
} else {
  tau_s = 25.0;
}

double k_s2 = G_k_s0 / k_k_s0; // eq. 22
double G_s = 1.0 * (1.0 - exp(-k_s2 * data->SNA)); // eq. 21 EDIT gain value changed to 1


// deltaHR_s
dY[51] = (-Y[51] + G_s * (data->SNA_buffer[(int)(T_s/DELTAT)-1]))/tau_s;

double tau_v;
if (data->PNA > data->PNA_buffer[0]){
  tau_v = 0.1;
} else {
  tau_v = 0.8;
}

double k_v2 = G_v0 / k_v0; // eq. 27
double G_v = 1.0 * (1.0 - exp(-k_v2 * data->PNA)); // eq. 26, EDIT gain value changed to 1 for more accurate description of HR used in

// deltaHR_v
dY[52] = (-Y[52] + G_v * (data->PNA_buffer[(int)(T_v/DELTAT)-1]))/tau_v;

// // SNA Effector Sites
// again the paper uses inconsistent notation
double T_sigma_lv = T_e_lv;
double T_sigma_rv = T_e_rv;
double T_sigma_V = T_e_V;
double T_sigma_R = T_e_R;

// The following are all described by eq. 29
// Left ventricular contractility
// sigma_lv
dY[53] = (-Y[53] + G_eff_lv * (data->SNA_buffer[(int)(T_sigma_lv/DELTAT)-1]))/tau_sigma_lv;

// Right ventricular contractility
// sigma_rv
dY[54] = (-Y[54] + G_eff_rv * (data->SNA_buffer[(int)(T_sigma_rv/DELTAT)-1]))/tau_sigma_rv;

// Venous tone
// sigma_V
dY[55] = (-Y[55] + G_eff_V * (data->SNA_buffer[(int)(T_sigma_V/DELTAT)-1]))/tau_sigma_V;

// Arterial resistance
// sigma_R
dY[56] = (-Y[56] + G_eff_R * (data->SNA_buffer[(int)(T_sigma_R/DELTAT)-1]))/tau_sigma_R;
