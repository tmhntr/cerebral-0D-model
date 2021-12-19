/*
Sanjay R. Kharche. Part of PM3 platforms.
Human circulation and dialysis unit model RHS. Sept 24, 2020.

Sept 24. 2020.
To extend this model using the baro-reflex model.
*/
static int f_ursino(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{

 realtype Y[NEQ];
 realtype dY[NEQ];
 int i;

 UserData data;
 data = (UserData) user_data;

  for(i=0;i<NEQ;i++) Y[i] = Ith(y,i+1);
  for(i=0;i<NEQ;i++) dY[i] = 0.0;
/*****************************************************************************************************************************/
/*****************************************************************************************************************************/
// %%% ALL PASSED PARAMETERS ARE LISTED BELOW %%%

double loc_t       		= data->p_ursino[1] - data->p_ursino[106]; // gnu fmod function. p11 is T_n_1, see above for definition of T_n_1.


data->p_ursino[105]   = loc_t; // passed back to main to calculate cardiac output.
double respRate  = data->p_ursino[12]; // % respiratory rate; units: (breaths/s)

#include "baroreflex.c"
double sigma_lv     				= Y[53];///1.5; //data->p_ursino[113]; // % units: mmHg/ml
double sigma_rv     				= Y[54];///0.935; //data->p_ursino[114]; //  units: mmHg/ml
// double delta_sigma_V 				= data->p_ursino[115]; //  units: ml
// double V_gain 							= 0.0; //dY[89]/2785.0; //delta_sigma_V / 2785; //  This is the change in venous volume / total venous volume
double sigma_R 							= Y[56]; //data->p_ursino[116]; //  units: mmHg/(ml/s)


/*****************************************************************************/
/**********************            FLOWS           **************************/
// % from Lim 2008, unless it is from Heldt 2002. Notation as close to the two papers as possible.
// % Flows into and out of each vascular compartment. Find units.
// % This structure works as a valve; allowing positive flow but not negative.

/************************ VENA CAVA ************************************/
// Y[18] is P_ra , Y[10] is P_sup
if(Y[10] > Y[18])	{
  data->Qsup    = (Y[10] - Y[18])/data->p_R[1];
}else {
  data->Qsup    = 0;
}

// Y[4] is P_ab , Y[11] is P_inf
data->Qab	      = (Y[4] - Y[11])/data->p_R[2]	;

// Y[11] is P_inf , Y[18] is P_ra
if(Y[11] > Y[18])	{
  data->Qinf	  = (Y[11] - Y[18])/data->p_R[3];
}	else {
  data->Qinf	  = 0;
}

/************************ R.A. OULET into R.V ******************************
Qrao is the flow coming out of the right atrium,
Rav is the atrio-ventricular value resistance. */
// Y[18] is P_ra , Y[12] is P_rv
if(Y[18] > Y[12]){
  data->Qrao = (Y[18] - Y[12])/data->p_R[4];
}	else {
  data->Qrao = 0;
}

/************************ R.V. OULET ******************************
Qro is the flow coming out of the right ventricle,
Rro is the pulmonary valve resistance.*/
// Y[12] is P_rv , Y[12] is P_pa
if (Y[12] > Y[13])	{
  data->Qro	 = (Y[12] - Y[13])	/data->p_R[5]	;
}	else {
  data->Qro	 = 0;
}

/********************  PULMONARY ARTERY ***************************************
Qpa is the flow through the pulmonary artery,
data->p_R[6] is the pulmonary artery resistance.*/
// Y[13] is P_pa , Y[14] is P_pv
data->Qpa	= (Y[13] - Y[14])	/data->p_R[6];

/********************* L.A. INLET / PULMONARY VEIN*****************************
Qli is the pulmonary venous flow AKA the flow into the L. Atrium,
data->p_R[7] is the pulmonary vein resistance. */
// Y[14] is P_pv , Y[17] is P_ra
if (Y[14] > Y[17]){
  data->Qli	 = (Y[14] - Y[17])	/data->p_R[7] ;
}else{
  data->Qli	= 0;
}

/************************ L.A. OULET -> L.V ***********************************
Qlao is the flow coming out of the left atrium,
Rmv is the mitral value resistance. */
// Y[17] is P_ra , Y[8] is P_rv
if(Y[17] > Y[8]) {
	data->Qlao = (Y[17] - Y[8]) / data->p_R[8];
}else {
  data->Qlao = 0.0;
}

/************************ L.V. OULET / AORTIC *********************************
Qlo is the flow out of the left ventricle,
Rlo is the aortic valve resistance. */
// Y[8] is P_rv , Y[9] is P_a
if (Y[8] > Y[9]){
  data->Qlo = (Y[8] - Y[9])/data->p_R[9];
}else	{
  data->Qlo		= 0;
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
// Y[9] is P_a , Y[0] is P_up
data->Qupi	= (Y[9] - Y[0])	/(data->p_R[13] * sigma_R)	;

// Y[0] is P_up, Y[10] is P_sup
if (Y[0] > Y[10]){
	data->Qupo	 = (Y[0]   - Y[10])	/data->p_R[14]	;
}else{
	data->Qupo= 0;
}

/*************************** SPLANCHNIC  ***********************************
Qsp1 is the splanchnic microcirculation?
Qsp2 is the splanchnic venous circulation,
Rsp1 & Rsp2 are the resistances of each.
Y[9] is P_a,
Y[2] is P_sp,
Y[4] is P_ab
*/
data->Qsp1 = (Y[9] - Y[2])	/(data->p_R[15] * sigma_R)	;

if(Y[4] < Y[2]){
  data->Qsp2 = (Y[2]   - Y[4])	/data->p_R[16]	;
}else{
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

data->Qk1 = (Y[9] - Y[1])/data->p_R[17];
data->Qk2 = (Y[1] - Y[4])/data->p_R[18];

/*************************** LOWER BODY ***********************************
Qll1 & Qll2 are the lower body/legs flow,
data->p_R[43] & data->p_R[44] are their resistances.
Y[12] is P_a,
Y[6] is P_ll,
Y[7] is P_ab.  */
data->Qll1	 	= (Y[9] - Y[3])	/(data->p_R[43]* sigma_R);

if (Y[3] > Y[4]){
	data->Qll2	= (Y[3]   - Y[4])	/data->p_R[44]	;
}else{
	data->Qll2	= 0;
}

/*******************************************************************************
Pumping function for the heart
*******************************************************************************/
#include "pumpingfxn.c"

#include "cerebral.c"


// // % eq. 12 from Ursino 2000 p206 (unless stated otherwise)
// // % mass transfer rate from ic compartment to is compartment
// double phi_k       = - k_K * ( ckic - beta_K * ckis ) ; // % Ursino 2008, eq 26. amount of solute exchanged at cellular membrane per unit time. units:
// double phi_na      = - k_Na * ( cnaic - beta_Na * cnais ) ; // % Ursino 2008, eq 26. amount of solute exchanged at cellular membrane per unit time. units:
// double phi_u       = - k_U * ( cuic - beta_U * cuis ) ; // % Ursino 2008, eq 26. amount of solute exchanged at cellular membrane per unit time. units:
// double phi_hco3    = - eta_hco3 * ( chco3ic - g_hco3 * chco3is ) ; // % Ursino 2000, eq 12. amount of solute exchanged at cellular membrane per unit time. units:
// double phi_h       = - eta_h * ( chic - g_h * chis ) ; // % Ursino 2000, eq 12. amount of solute exchanged at cellular membrane per unit time. units:
// double phi_p       = 0.0; // % Ursino 2000. p207, second para.
// double phi_cl      = phi_na + phi_k + phi_h - phi_hco3; // % Ursino 2000, eq 13.
// // %
// // % eq 11 from Ursino 2000. These are eq 20 to 22 p207 of Ursino 2000. Algebriac equations.
// double Rhco3ic 	= etaprime_r * ( kprime_a * cco2ic - chco3ic * chic ); // % eq 20
// // % for chpic in eq 21, see bottom para p207 of Ursino 2000.
// double chpic    = cpic0 - cpic; // % ref ursino 2000 p. 207 last paragraph
// double Rpic 		= etaprimeprime_r * (kprimeprime_a * chpic - cpic * chic ); // % eq 21. Dont know chpic yet (4 June 2020).
// double Rhic 		= Rhco3ic + Rpic; // % eq 22.
// // %
// // % Rsex values are analogous to Rsic values. Ursino 2000 p.
// double Rhco3is 	= etaprime_r * ( kprime_a * cco2ex - chco3is * chis ); // % eq 20
// // %
// // % for chpis in eq 21, see bottom para p207 of Ursino 2000.
// // % Protein activity is assumed to be mainly in plasma
// double chpis    = c_pis - cpis; // % ref ursino 2000 p. 207 last paragraph
// double Rpis 		= etaprimeprime_r * (kprimeprime_a * chpis - cpis * chis ); // % etaprimeprime_r * (kprimeprime_a * chpis - cpex * chex ); // % eq 21. Dont know chpic yet (4 June 2020).
// double Rhis 		= Rhco3is + Rpis; // % eq 22.
// double Rhco3pl 	= etaprime_r * ( kprime_a * cco2ex - chco3pl * chpl ); // % eq 20
// // %
// // % for chpic in eq 21, see bottom para p207 of Ursino 2000.
// double chppl    = c_ppl - cppl; // % ref ursino 2000 p. 207 last paragraph
// double Rppl 		= etaprimeprime_r * (kprimeprime_a * chppl - cppl * chpl ); // % eq 21. Dont know chpic yet (4 June 2020).
// double Rhpl 		= Rhco3pl + Rppl; // % eq 22.
// // %
// // % TEMPORARY reaction rate was entered as the sum of rates in pl and is
// // % compartments. Find reference or revise
// double Rhco3ex  = Rhco3pl + Rhco3is; double Rpex = Rppl + Rpis; double Rhex = Rhpl + Rhis;
// // %
// double V_rc     = 1.3;          // % red blood cell volume, units (L) ursino + innocenti 2008
// double V        = V_rc + Y[3-1];  // % whole blood volume
// double Hct      = V_rc / V;
//
// // % eq 29 from Ursino and Innocenti 2008 p888
// double Q_eK 		 = Q_B*(F_p * (1 - Hct) + F_R * gamma_K 	* alphac);
// double Q_eNa	   = Q_B*(F_p * (1 - Hct) + F_R * gamma_Na 	* alphac);
// double Q_eU		   = Q_B*(F_p * (1 - Hct) + F_R * gamma_U 	* R_DU);
// double Q_eHCO3	 = Q_B*(F_p * (1 - Hct) + F_R * gamma_HCO3 * alphac);
// double Q_eCl		 = Q_B*(F_p * (1 - Hct) + F_R * gamma_Cl 	* alphaa);
//
// // % eq 28 from Ursino and Innocenti 2008 p888
// // % convective and diffusive transport to dialyzer.
// // J is solute removal rate across the dialyser.
// double J_k     	= (D_s * 		(1.0 - QF / Q_eK) + QF) 			* ckex 	- D_s * 	     	(1.0 - QF / Q_eK) 	* ckd ;
// double J_na    	= (D_s * 		(1.0 - QF / Q_eNa) + QF) 		* cnaex 	- D_s *  		(1.0 - QF / Q_eNa) 	* cnad ;
// double J_u     	= (D_u * 		(1.0 - QF / Q_eU) + QF) 		* cuex 	- D_u *  		(1.0 - QF / Q_eU) 	* cud ;
// double J_hco3  	= (D_hco3 * 	(1.0 - QF / Q_eHCO3) 	+ QF) 	* chco3ex- D_hco3 *     	(1.0 - QF / Q_eHCO3) * chco3d ;
// double J_cl    	= (D_s * 		(1.0 - QF / Q_eCl) 		+ QF) 	* cclex 	- D_s * 	      	(1.0 - QF / Q_eCl) 	* ccld ;
// double J_h      = 0.0; // % for hydrogen J is 0.
// double J_p      = 0.0; // % for proteins, J is 0.


/***********************ODES***************************************************/
/******************************************************************************/
// %================== ALL ODEs ARE LISTED BELOW ==================

// % in Urisino eq 2, Oic is the intracellular osmotic activity. units: probably mEq/L or mmol/L unless the 0.93 has
// % Ois is the interstitial osmotic activity.
// dY[0] 	=   k_f*(cic - cis); // % equation 1, appendix 1, Vic, intracellular fluid volume. Ursino 2000 units: L (liters).
// // is the 1000 conversion factor applied consistently throughout the RHS?
// dY[1] 	= - k_f*(cic - cis) + (Fa - Rv)/1000.; // % eq. 2, appendix 1, Vis. Vis is interstial fluid volume. Ursino 2000. units: L.
// dY[2] 	= ( -(Fa - Rv) - QF + Qinfused)/1000.; 	// % eq. 3, appedix 1, Vpl. Vpl is plasma volume. Ursino 2000. unit: L.

/* Y[3] is Pup, pressure in upper body, head and neck.
 		eq. 22, appendix 2, p585.*/
dY[0]     = (data->Qupi - data->Qupo) / data->p_C[10];


/*************** kidney ************************************/
/*************** original formula **************************/
dY[1]     = (data->Qk1 - data->Qk2) / data->p_C[13]; // % Y[5-1] is Pk, Ck is non-zero. eq. 23, appendix 2, p585.

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
dY[2]     = (data->Qsp1 - data->Qsp2) / data->p_C[11];

// % Y[6] is Pll, legs pressure.
dY[3]     = (data->Qll1 - data->Qll2) / data->p_C[26];

// % Y[7] is Pab.
dY[4]     = (data->Qk2 + data->Qsp2 + data->Qll2 - data->Qab) / data->p_C[2] ;

// % Y[8] is Pth: thoracic pressure. aka Y8 is Pbias. p1241, Heldt 2002 paper.
dY[5]     = (2.0*M_PI*respRate)*cos(2.0*M_PI*respRate*data->p_ursino[1]); // % This varies with respiratory variation, mean value taken from Heldt thesis 2.3 p.47 resprate units are (breath/s).

/***************************** ODEs for Compliances **************************/
dY[6] = (Clv - Y[6])/DELTAT;    // Left Ventricle Compliance
dY[7] = (Crv - Y[7])/DELTAT;  // Right Ventricle Compliance
dY[15] = (Cla - Y[15])/DELTAT;  // Left Atrial Compliance
dY[16] = (Cra - Y[16])/DELTAT;  // Right Atrial Compliance

/****************************************************************************************************************************/
// % Y[11] is LV pressure.
dY[8]    = (1.0 / Y[6])  * ((Y[5] - Y[8]) * dY[6] + data->Qlao - data->Qlo ) + dY[5];

// Y34 is left atrial pressure.
dY[17]    = (1.0 / Y[15]) * ((Y[5] - Y[17]) * dY[15] + data->Qli - data->Qlao) + dY[5];

// Y[15] is Right Ventricle Pressure
dY[12]    = (1.0 / Y[7]) * ((Y[5] - Y[12]) * dY[7] + data->Qrao - data->Qro) + dY[5];

// Y[35] is right atrial pressure
dY[18]    = (1.0 / Y[16]) * ((Y[5] - Y[18]) * dY[16] + data->Qinf + data->Qsup - data->Qrao) + dY[5];

// % Y[12] is Pa, ascending aorta pressure.
// if(data->argv[3] == 1.0){
//   dY[12]    = (1.0/data->p_C[6]) * (data->Qlo + ( -(Fa - Rv) - QF + Qinfused) - (data->Q_ao_bra + data->Q_ao_tho)) + 0.333*dY[8];
// }else if(data->argv[3] == 0.0){
//   dY[12]    = (1.0/data->p_C[6]) * (data->Qlo - (data->Q_ao_bra + data->Q_ao_tho)) + 0.333*dY[8];
// }
// Y[12] is systemic arterial pressure (aorta)
dY[9]    = (1.0/data->p_C[6]) * (data->Qlo - (data->Qupi + data->Qsp1 + data->Qk1 + data->Qll1 + Q_Ceri)) + 0.333*dY[5];

// % Y[13] is Psup, eq. 27 appendix 2.
dY[10]    = ((data->Qupo + Q_Cero) - data->Qsup) / data->p_C[1] + dY[5];
// % Y[14] is Pinf, eq. 28 appendix 2.
dY[11]    = (data->Qab - data->Qinf) / data->p_C[3] + dY[5];

// % Y[16] is Ppa. Pulmonary circulation inlet pressure.
dY[13]    = (data->Qro - data->Qpa) / data->p_C[4] + dY[5]; // % Y[17-1] is Ppa, eq. 30 appendix 2.

// % Y[17] is Ppv. Left ventricle inlet pressure.
dY[14]    = (data->Qpa - data->Qli) / data->p_C[5] + dY[5]; // % Y[18-1] is Ppv, eq. 31 appendix 2.

// //  Brachiocephalic Aortic Pressure
// dY[50] = (data->Q_ao_bra - (data->Qupi + Q_Ceri))/data->p_C[7] + dY[8]; //  Brachiocephalic Aortic Pressure
//
// //  Thoracic Aortic Pressure
// dY[51] = (data->Q_ao_tho - data->Q_ao_abd)/data->p_C[8] + dY[8];
//
// //  Abdominal Aortic Pressure
// dY[52] = (data->Q_ao_abd - (data->Qsp1 + data->Qll1 + data->Qk1))/data->p_C[9] + dY[8];


/****************************************************************************************************************************/
// % for eq 11, Ursino 2000
// % State variables 19-25 are Intracellular Solute mass
// dY[19-1] = phi_u;                   // % urea 		   eq 26 Ursino 2000.
// dY[20-1] = phi_na;                  // % sodium 	   eq 26 Ursino 2000.
// dY[21-1] = phi_k;                   // % potassium   eq 26 Ursino 2000.
// dY[22-1] = phi_cl;                  // % chlorine 	 eq 11 Ursino 2000.
// dY[23-1] = phi_hco3 + Rhco3ic;      // % bicarbonate eq 11 Ursino 2000.
// dY[24-1] = phi_h 		+ Rhic;         // % hydrogen 	 eq 11 Ursino 2000.
// dY[25-1] = phi_p 		+ Rpic;         // % urea 		   eq 11 Ursino 2000.
//
// // % for eq 14, Ursino 2000
// // % State variables 28-34 are Extracellular Solute mass
// dY[26-1]   = -phi_u 	  - J_u/1000. 		+ (Qinfused/1000.)*cuinf; // % urea eq 14, p207, Ursino 2000. combined with eqn 27, Ursino 2008.
// dY[27-1]   = -phi_na    - J_na/1000.		+ (Qinfused/1000.)*cnainf; // % sodium 	eq 14, p207, Ursino 2000.combined with eqn 27, Ursino 2008.
// dY[28-1]   = -phi_k 	  - J_k/1000. 		+ (Qinfused/1000.)*ckinf; // % potassium  eq 14, p207, Ursino 2000.combined with eqn 27, Ursino 2008.
// dY[29-1]   = -phi_cl 	  - J_cl/1000.    + (Qinfused/1000.)*cclinf; // % chlorine 	eq 14, p207, Ursino 2000.
// dY[30-1]   = -phi_hco3 	- J_hco3/1000. 	+ Rhco3ex; // % bicarbonate eq 14, p207, Ursino 2000.
// dY[31-1]   = -phi_h 	  - J_h/1000. 		+ Rhex	; // % hydrogen eq 14, p207, Ursino 2000.
// dY[32-1]   = -phi_p 	  - J_p/1000. 		+ Rpex	; // % protein eq 14, p207, Ursino 2000.
/*****************************************************************************************************************************/
/******************************************************************************************************************************/


// printf("%d=%f\n",57, Y[56]);
//
// printf("C_dml=%f\n", data->C_dqs[0]);
// printf("dC_dml=%f\n", dC_dml);
// printf("q_dml=%f\n", data->q_ml);


 for(i=0;i<NEQ;i++) Ith(ydot,i+1) = dY[i];
 return(0);
}
