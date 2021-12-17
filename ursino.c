/* Sanjay R. Kharche.
The combined Heldt-dialysis model, nominally called Ursino model. This model combines
Heldt whole body circulation and Ursino dialysis unit into a 0D discription.
The names of these files are misnomers, if anything, they should be called Lim and Eu Bo Shim.
References:
The composite model is desciribed in this paper. In the pdfs/ directory, refs refer to references in this paper:
1) Ki Moo Lim, Sung Wook Choi, Byung Goo Min, and Eun Bo Shim.  Numerical Simulation of the Effect of Sodium Profile on Cardiovascular Response to Hemodialysis. Yonsei Med J 49(4):581 - 591, 2008.
2) Heldt 2002. Journal of Applied Physiology.
3) shear from Secomb.
4) detailed kidney as a generalization of single resisitance-capacitor.

Further, papers by Secomb have been used to implement myogenic response and shear.

Units in this model are:
time: seconds.
volume: Liters.
flow: ml/s.
pressure: mmHg.
concentration: mmol
resistance: mmHg*s/ml
capacitance: ml/mmHg
elastance:
*/
// my standard headers and functions.
#include "ursino.h"

int main (int argc, char *argv[])
{
// exit(-1);
void *cvode_mem;  // pointer to memory: the full state lives here.
realtype t, tout;
int iout, NOUT, retval, i;
char *str;
FILE *cerebral,*cardiac_output_file,*pinkFile,*expFile;
UserData 	data; // instance pointer.
data 		= (UserData) malloc(sizeof *data); // now it is created. // allocated memory to pointer.

/* Create serial vector of length NEQ for I.C. and abstol */
N_Vector 	y_ursino = N_VNew_Serial(NEQ); // allocated memory to pointer.

// Ith(y_ursino, 1)   = 25;    		// % Vic, 	units: L
// Ith(y_ursino, 2)   = 11;    		// % Vis, 	units: L
// Ith(y_ursino, 3)   = 3.25;  		// % Vpl, 	units: L
Ith(y_ursino, 0 + 1)   = 13.7;  		// % Pup, 	units: mmHg
Ith(y_ursino, 1 + 1)   = 8.0;  		// % Pk, 		units: mmHg. This is now 12 variables.
Ith(y_ursino, 2 + 1)   = 8.3;     // % Psp, 	units: mmHg
Ith(y_ursino, 3 + 1)   = 8.3;  		// % Pll, 	units: mmHg
Ith(y_ursino, 4 + 1)   = 6.2;  		// % Pab, 	units: mmHg
Ith(y_ursino, 5 + 1)   = -5.0;    		// % Pth, 	units: mmHg. // % is pth a variable, or a constant? 5 Sept. 2020.
Ith(y_ursino, 6 + 1)  = 10.0;  		// % Cl, 	units: ml/mmHg
Ith(y_ursino, 7 + 1)  = 20.0;  		// % Cr, 		units: ml/mmHg
Ith(y_ursino, 8 + 1)  = 12.78;  		// % Pl, left ventricular pressure, units: mmHg
Ith(y_ursino, 9 + 1)  = 80.0;  			// % Pa, aortic pressure, units: mmHg.
Ith(y_ursino, 10 + 1)  = 5.0;     	// % Psup, 		units: mmHg
Ith(y_ursino, 11 + 1)  = 5.0;     	// % Pinf, 		units: mmHg
Ith(y_ursino, 12 + 1)  = 5.60;   		// % Pr, 			units: mmHg. right ventricle pressure.
Ith(y_ursino, 13 + 1)  = 25.66;  		// % Ppa,			units: mmHg. Pulmonary artery pressure.
Ith(y_ursino, 14 + 1)  = 10.99;  		// % Ppv, 		units: mmHg. Pulmonary vein pressure.
// Ith(y_ursino, 19)  = 100.0;   		// % Muic, 		units: mmol
// Ith(y_ursino, 20)  = 250.0; 		// % Mnaic, 	units: mmol
// Ith(y_ursino, 21)  = 3535.0; 		// % Mkic, 		units: mmol
// Ith(y_ursino, 22)  = 84.0;  		// % Mclic, 	units: mmol
// Ith(y_ursino, 23)  = 10.0;  		// % MHco3ic, units: mmol
// Ith(y_ursino, 24)  = 100.0; 		// % Mhic,	 	units: mmol
// Ith(y_ursino, 25)  = 0.0;   		// % Mpic, 		units: mmol
// Ith(y_ursino, 26)  = 55.0;  		// % Muex, 		units: mmol
// Ith(y_ursino, 27)  = 2130.0; 		// % Mnaex, 	units: mmol
// Ith(y_ursino, 28)  = 75.0;  		// % Mkex, 		units: mmol
// Ith(y_ursino, 29)  = 1470.0; 		// % Mclex, 	units: mmol
// Ith(y_ursino, 30)  = 100.0; 		// % Mhco3ex, units: mmol
// Ith(y_ursino, 31)  = 100.0; 		// % Mhex, 		units: mmol
// Ith(y_ursino, 32)  = 0.0;   		// % Mpex, 		units: mmol

// l/r atria
Ith(y_ursino, 15 + 1)  = 1.0; 			// Cla variable, left atrial elastance.
Ith(y_ursino, 16 + 1)  = 2.0;   		//  Cra variable, right atrial elastance.

Ith(y_ursino, 17 + 1)  = 8.0;   		//  left atrial pressure initial condition, mmHg.
Ith(y_ursino, 18 + 1)  = 1.0;   		//  right atrial pressure initial condition, mmHg.


/**** detailed kidney not currently in use ********/
// all kidney pressures, including inlet right and left pressures.
// for(i=37;i<=48;i++) Ith(y_ursino, i) = 8.0;
// Ith(y_ursino, 49)  = 70.0;   //  Right Kidnet Inlet Pressure, mmHg.
// Ith(y_ursino, 50)  = 65.0;   //  Left Kidnet Inlet Pressure, mmHg.
/********************************************************/

//
// Ith(y_ursino, 51)  = 80.0;   //  Brachiocephalic Aortic Pressure
// Ith(y_ursino, 52)  = 80.0;   //  Thoracic Aortic Pressure
// Ith(y_ursino, 53)  = 92.0;   //  Abdominal Aortic Pressure

//*******************************************************************************************************************************
// Cerebral Pressures
// for(i=54;i<=65;i++) Ith(y_ursino, i) = 5.0;
Ith(y_ursino, 19 + 1)  = 9.5; // eq 1: dP_ic
Ith(y_ursino, 20 + 1)  = 25.0; // P_c
Ith(y_ursino, 21 + 1)  = 14.0; // eq. 11 : dP_v
Ith(y_ursino, 22 + 1)  = 58.5; // dP_djs
Ith(y_ursino, 23 + 1)  = 58.5; // dP_djs
Ith(y_ursino, 24 + 1)  = 58.5; // dP_djs
Ith(y_ursino, 25 + 1)  = 58.5; // dP_djs
Ith(y_ursino, 26 + 1)  = 58.5; // dP_djs
Ith(y_ursino, 27 + 1)  = 58.5; // dP_djs
Ith(y_ursino, 28 + 1)  = 92.5; // P_ICAl
Ith(y_ursino, 29 + 1)  = 92.5; // P_ICAr
Ith(y_ursino, 30 + 1)  = 92.5; // P_BA
//
// Ith(y_ursino, 66)  = 10.1 * 0.30; // V_djs
// Ith(y_ursino, 67)  = 10.1 * 0.12; // V_djs
// Ith(y_ursino, 68)  = 10.1 * 0.08; // V_djs
// Ith(y_ursino, 69)  = 10.1 * 0.30; // V_djs
// Ith(y_ursino, 70)  = 10.1 * 0.12; // V_djs
// Ith(y_ursino, 71)  = 10.1 * 0.08; // V_djs
// Ith(y_ursino, 72)  = 4.5; // V_vi

for(i=31;i<=42;i++) {
	Ith(y_ursino, i + 1) = 0.0; // X_aut_djs + X_
}

for(i=43;i<=48;i++) {
	Ith(y_ursino, i + 1) = 0.004; // C_djs
}

Ith(y_ursino, 49 + 1)  = 90.0; // P_aff
Ith(y_ursino, 50 + 1)  = 90.0; // X0 for P_error
Ith(y_ursino, 51 + 1)  = 0.0; // deltaHR_s
Ith(y_ursino, 52 + 1)  = 0.0; // deltaHR_v
Ith(y_ursino, 53 + 1)  = 1.0; // sigma_lv
Ith(y_ursino, 54 + 1)  = 1.0; // sigma_rv
Ith(y_ursino, 55 + 1)  = 0.0; // sigma_V
Ith(y_ursino, 56 + 1)  = 1.0; // sigma_R

//*******************************************************************************************************************************

cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
CVodeInit(cvode_mem, f_ursino, 0.0, y_ursino);
CVodeSStolerances(cvode_mem, ATOL, RTOL);
CVDense(cvode_mem, NEQ);
CVodeSetMaxStep(cvode_mem,DELTAT);

#include "p_ursino.c"
// if(0==0) printf("%d : My number of Arguments = %d,\n", atoi(argv[1]), argc);

//*******************************************************************************************************************************
//**************************** END OF BAROREFLEX PARAMETERS *********************************************************************
//*******************************************************************************************************************************


// opening various files for writing
// In future can combine multiple files to save computation time

str = malloc(32*sizeof(char)); sprintf(str,"cerebral_%d_%d_%d.dat", atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
cerebral = fopen(str,"w+"); free(str);
//
// str = malloc(32*sizeof(char)); sprintf(str,"flows%05d.dat", atoi(argv[1]));
// flows = fopen(str,"w+"); free(str);
//
str = malloc(32*sizeof(char)); sprintf(str,"pressures_%d_%d_%d.dat", atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
cardiac_output_file = fopen(str,"w+"); free(str);

int headersPrinted = 0;

// printf("Reading noise files...");
str = malloc(128*sizeof(char)); sprintf(str,"samples/pinkNoise_%05d.dat", atoi(argv[1]));
pinkFile = fopen(str, "r");
free(str);

double pinkNoise[numBeats] = {0.0};
double expRand[numBeats] = {0.0};

for (int i = 0; i < numBeats; i++){
    fscanf(pinkFile, "%lf", &pinkNoise[i]);
}
fclose(pinkFile);

if (atoi(argv[2]) == 1){
	str = malloc(128*sizeof(char)); sprintf(str,"samples/expRand_%05d.dat", atoi(argv[1]));
	expFile = fopen(str, "r");
	free(str);


	for (int i = 0; i < numBeats; i++){
	    fscanf(expFile, "%lf", &expRand[i]);
	}
	fclose(expFile);
}

// printf("Read complete.");


// JJ TH Nov 3
// calculate length of each buffer
int s_l = (int)(data->p_ursino[151]/DELTAT);
int p_l = (int)(data->p_ursino[143]/DELTAT);

// Declaring Buffers
// double SNA_buffer[s_l], PNA_buffer[p_l];
// for (int s = 0; s < s_l; s++){
// 	SNA_buffer[s] = 3.0;
// }
// for (int p = 0; p < p_l; p++){
// 	PNA_buffer[p] = 3.0;
// }
double deltaHR;


//***********************************************************************
// initialise iterators.

// NOUT = (int)(SIMTIME/DELTAT);


iout = 0; tout = DELTAT;
int cardiac_iter = 0;
double cardiac_output = 0.0; // this is time integral of Qlo (LV output flow) over one heart period. for now, it is p11 as the reflex is switched off 10 Oct 2020.

double flowpermin[7]		 		= {0.0};

double sysPressures[7] 			= {0.0};
double PresToCompare[7] 		= {0.0};
double diasPressures[7] 		= {0.0};
double maxFlows[7] 					= {0.0};
double FlowsToCompare[7] 		= {0.0};
double minFlows[7] 					= {0.0};

// printf("Starting solver loop...");
// TIME LOOP STARTS HERE
while(cardiac_iter<numBeats){

  data->p_ursino[1] = tout; // pursino1 is time in seconds. time dependent paramter.

  // ***************************************************************************************************
  // Baroreflex Implementation JJ TH Nov 3
  // #include "BR_Integration.c"

  // Update SNA_buffer
  for(int s = s_l-1; s>0; s--) data->SNA_buffer[s] = data->SNA_buffer[s-1];
  data->SNA_buffer[0] = data->SNA;

  // Update PNA_buffer
  for(int p = p_l-1; p>0; p--) data->PNA_buffer[p] = data->PNA_buffer[p-1];
  data->PNA_buffer[0] = data->PNA;

  // #include "BR_Effector.c"


  // data->p_ursino[11] = (60.0/HR);

  // printf("heart rate calculations");
  // Timing Parameters
  if(data->p_ursino[1] >= data->p_ursino[106] + data->RR[0]){

  	data->p_ursino[106]  		= data->p_ursino[1];
		data->RR[2] 						= data->RR[1];
		data->RR[1] 						= data->RR[0];

		deltaHR = (19.64*Ith(y_ursino, 51 + 1)) - (17.95*Ith(y_ursino, 52 + 1)) - (1.225*pow(Ith(y_ursino, 51 + 1),2)) + (1.357*pow(Ith(y_ursino, 52 + 1),2)) - (1.523*Ith(y_ursino, 51 + 1)*Ith(y_ursino, 52 + 1));
		data->RR[0] = 60.0/(data->p_ursino[148+19]/* + deltaHR*/);

		if (atoi(argv[2]) < 1){
			data->RR[0] 				= data->RR[0] + pinkNoise[cardiac_iter]*(data->p_ursino[169]*(data->RR[0])/0.033);
		} else {
			// The following is AF condition 1 and 3 as described by Scarsoglio et al. 2014
			data->p_ursino[171] = expRand[cardiac_iter]/(-9.2*data->RR[0] + 14.6); //r4_exponential_sample(1.0/(-9.2*data->RR[0] + 14.6));
			data->p_ursino[170] = (data->RR[0] - 1.0/(-9.2*data->RR[0] + 14.6)) + pinkNoise[cardiac_iter]*(data->p_ursino[169]*(data->RR[0] - 1.0/(-9.2*data->RR[0] + 14.6))/0.033);

			data->RR[0] 				= data->p_ursino[171] + data->p_ursino[170];

			// printf("RR: %f\ntau: %f\nnu: %f\n\n",data->RR[0], data->p_ursino[170], data->p_ursino[171]);

			// EQ 2 from Scarsoglio et al. 2014
			data->Esys_lv 			= 0.59 * (data->RR[1]/data->RR[2]) + 0.91;
		}
		// HR = data->p_ursino[148+19];// + deltaHR; // eq. 18 of paper ABC.

  	data->Ts             		= 0.37	*sqrt(data->RR[0]);
  	data->Tasys          		= 0.25	*sqrt(data->RR[0]);
  	data->Tav            		= 0.19	*sqrt(data->RR[0]);
		cardiac_iter++;
  }
  // END OF BAROREFLEX OPERATIONS
  // ***************************************************************************************************


	// CEREBRAL AUTOREGULATION
	// ***************************************************************************************************

	// #include "autoregulation.c"
	// printf("q_dml=%f\t X_aut=%f\t X_CO2=%f\t A=%f\t C_ml=%f\n", data->q_ml, data->x_aut[1], data->x_CO2[1], data->A_CO2[1], Ith(y_ursino, 66));
	//
	// printf("q_dal=%f\t X_aut=%f\t X_CO2=%f\t A=%f\t C_al=%f\n", data->q_al, data->x_aut[2], data->x_CO2[2], data->A_CO2[2], Ith(y_ursino, 67));


	// END OF CEREBRAL AUTOREGULATION
	// ***************************************************************************************************

  CVodeSetUserData(cvode_mem, data); // you must tell cvode_mem about data. You have time dependent data.
  retval = CVode(cvode_mem, tout, y_ursino, &t, CV_NORMAL); if(check_retval(&retval, "CVode", 1)) exit(1);

// JJJ Shear Stress

// array of lengths for each flow. Might be better to choose reasonable lengths. But for now all set to Sanjay's constant parameter.
for (i = 0; i<45; i++){
  data->lengths[i] = data->p_ursino[101];
}

// mu is p[100].
for(i = 1; i<45; i++){
  data->radius[i]    	= pow(( (8.0 * data->lengths[i]  * data->p_ursino[100] )/ (M_PI * data->p_R[i]) ), 0.25 );
  data->tau[i]       	= 4.0 * data->p_ursino[100] * (*(data->flows[i]))  / (M_PI * pow(data->radius[i],3));
}

  // Configure output files:
  #include "output_ursino.c"

  tout = tout + DELTAT;
} // end of time loop.


	// for (i = 1; i < 72; i++){
	// 	fprintf(cardiac_output_file, "%f\t", Ith(y_ursino, i));
	// }
	// for (i = 0; i < 6; i++){
	// 	fprintf(cardiac_output_file, "%f\t%f\t%f\t%f\t%f\t", data->x0_aut[i], data->x_aut[i], data->A_CO2[i], data->x0_CO2[i], data->x_CO2[i]);
	// }
	fclose(cerebral);
	// fclose(flows);
  fclose(cardiac_output_file);

  N_VDestroy_Serial(y_ursino);
	// emxDestroyArray_real_T(X);
  CVodeFree(&cvode_mem);
  free(data);
	printf("done run %d\n",atoi(argv[1]));

	return 0; // the main must return a success exit code to the middleware.
}
