/* Sanjay R. Kharche.
The combined Heldt-dialysis model, nominally called Ursino model. This model combines
Heldt whole body circulation and Ursino dialysis unit into a 0D discription.
The names of these files are misnomers, if anything, they should be called Lim and Eu Bo Shim.
References:
The composite model is desciribed in this paper. In the pdfs/ directory, refs refer to references in this paper:
1) Ki Moo Lim, Sung Wook Choi, Byung Goo Min, and Eun Bo Shim.  Numerical Simulation of the Effect of Sodium Profile on Cardiovascular Response to Hemodialysis. Yonsei Med J 49(4):581 - 591, 2008.
2) Heldt 2002. Journal of Applied Physiology.

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

char *str, *stateFilename, *outputInfoFilename, *postprocessedFilename;
FILE *stateFile, *outputInfoFile, *postprocessedFile, *pinkFile, *expFile;
UserData 	data; // instance pointer.
data 		= (UserData) malloc(sizeof *data); // now it is created. // allocated memory to pointer.

/* Create serial vector of length NEQ for I.C. and abstol */
N_Vector 	y_ursino = N_VNew_Serial(NEQ); // allocated memory to pointer.



//******************************************************************************
//*** Heldt model **************************************************************
//******************************************************************************

// Blood pressure states
Ith(y_ursino, 0 + 1)   = 20.7;  		// % Pup, 	units: mmHg
Ith(y_ursino, 1 + 1)   = 20.0;  		// % Pk, 		units: mmHg.=
Ith(y_ursino, 2 + 1)   = 20.0;     // % Psp, 	units: mmHg
Ith(y_ursino, 3 + 1)   = 20.0;  		// % Pll, 	units: mmHg
Ith(y_ursino, 4 + 1)   = 5.2;  		// % Pab, 	units: mmHg

// Thoracic pressure (pressure in thoracic cavity NOT blood pressure)
Ith(y_ursino, 5 + 1)   = -5.0;    		// % Pth, 	units: mmHg.

// Ventricular states (BP and compliances)
Ith(y_ursino, 6 + 1)  = 10.0;  		  // % Cl, 	units: ml/mmHg
Ith(y_ursino, 7 + 1)  = 20.0;  		  // % Cr, 		units: ml/mmHg
Ith(y_ursino, 8 + 1)  = 12.78;  		// % Pl, left ventricular pressure, units: mmHg

// More blood pressure states
Ith(y_ursino, 9 + 1)  = 80.0;  		// % Pa, aortic pressure, units: mmHg.
Ith(y_ursino, 10 + 1)  = 5.0;     	// % Psup, 		units: mmHg
Ith(y_ursino, 11 + 1)  = 5.0;     	// % Pinf, 		units: mmHg
Ith(y_ursino, 12 + 1)  = 5.60;   		// % Pr, 			units: mmHg. right ventricle pressure.
Ith(y_ursino, 13 + 1)  = 25.66;  		// % Ppa,			units: mmHg. Pulmonary artery pressure.
Ith(y_ursino, 14 + 1)  = 10.99;  		// % Ppv, 		units: mmHg. Pulmonary vein pressure.

// Atrial states (BP and compliances)
Ith(y_ursino, 15 + 1)  = 1.0; 			// Cla variable, left atrial elastance.
Ith(y_ursino, 16 + 1)  = 2.0;   		//  Cra variable, right atrial elastance.
Ith(y_ursino, 17 + 1)  = 8.0;   		//  left atrial pressure initial condition, mmHg.
Ith(y_ursino, 18 + 1)  = 1.0;   		//  right atrial pressure initial condition, mmHg.
//******************************************************************************



//******************************************************************************
//*** Ursino and Giannessi (2010) model ****************************************
//******************************************************************************

// Cerebral Pressures
Ith(y_ursino, 19 + 1)  = 9.5; // eq 1: dP_ic
Ith(y_ursino, 20 + 1)  = 25.0; // P_c
Ith(y_ursino, 21 + 1)  = 14.0; // eq. 11 : dP_v
Ith(y_ursino, 22 + 1)  = 50.5; // P_djs
Ith(y_ursino, 23 + 1)  = 50.5; // P_djs
Ith(y_ursino, 24 + 1)  = 50.5; // P_djs
Ith(y_ursino, 25 + 1)  = 50.5; // P_djs
Ith(y_ursino, 26 + 1)  = 50.5; // P_djs
Ith(y_ursino, 27 + 1)  = 50.5; // P_djs
Ith(y_ursino, 28 + 1)  = 85.5; // P_ICAl
Ith(y_ursino, 29 + 1)  = 85.5; // P_ICAr
Ith(y_ursino, 30 + 1)  = 85.5; // P_BA

// Cerebral autoregulation states
for(i=31;i<=42;i++) {
  Ith(y_ursino, i + 1) = 0.0; // X_aut_djs, X_co2_djs
}

// Cerebral compliance states
for(i=43;i<=48;i++) {
  Ith(y_ursino, i + 1) = 0.004; // C_djs
}
//******************************************************************************



//******************************************************************************
//*** Lin model ****************************************************************
//******************************************************************************

// Baroreflex states
Ith(y_ursino, 49 + 1)  = 90.0; // P_aff
Ith(y_ursino, 50 + 1)  = 90.0; // X0 for P_error
Ith(y_ursino, 51 + 1)  = 0.0; // deltaHR_s
Ith(y_ursino, 52 + 1)  = 0.0; // deltaHR_v
Ith(y_ursino, 53 + 1)  = 1.0; // sigma_lv
Ith(y_ursino, 54 + 1)  = 1.0; // sigma_rv
Ith(y_ursino, 55 + 1)  = 0.0; // sigma_V
Ith(y_ursino, 56 + 1)  = 1.0; // sigma_R

//******************************************************************************



//******************************************************************************
//*** Solver setup *************************************************************
//******************************************************************************

cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
CVodeInit(cvode_mem, f_ursino, 0.0, y_ursino);
CVodeSStolerances(cvode_mem, ATOL, RTOL);
CVDense(cvode_mem, NEQ);
CVodeSetMaxStep(cvode_mem,DELTAT);

#include "p_ursino.c"

//******************************************************************************



//******************************************************************************
//*** Input/output setup *******************************************************
//******************************************************************************

stateFilename = malloc(32*sizeof(char));
sprintf(stateFilename,"statesOutput.%d.dat", atoi(argv[1])); // argv2
stateFile = fopen(stateFilename,"w+");

outputInfoFilename = malloc(32*sizeof(char));
sprintf(outputInfoFilename,"outputInfo.%d.txt", atoi(argv[1])); // argv2
outputInfoFile = fopen(outputInfoFilename,"w+");

postprocessedFilename = malloc(32*sizeof(char));
sprintf(postprocessedFilename,"result.%d.txt", atoi(argv[1])); // argv2
postprocessedFile = fopen(postprocessedFilename,"w+");



free(stateFilename);
free(outputInfoFilename);
free(postprocessedFilename);

int headersPrinted = 0;
//******************************************************************************



//******************************************************************************
//*** Parse input files ********************************************************
//******************************************************************************

str = malloc(128*sizeof(char)); sprintf(str,"samples/pinkNoise_%05d.dat", atoi(argv[1]));
pinkFile = fopen(str, "r");
free(str);

double pinkNoise[numBeats] = {0.0};

for (int i = 0; i < numBeats; i++){
    fscanf(pinkFile, "%lf", &pinkNoise[i]);
}
fclose(pinkFile);

double expRand[numBeats] = {0.0};
if (atoi(argv[2]) == 1){
	str = malloc(128*sizeof(char)); sprintf(str,"samples/expRand_%05d.dat", atoi(argv[1]));
	expFile = fopen(str, "r");
	free(str);


	for (int i = 0; i < numBeats; i++){
	    fscanf(expFile, "%lf", &expRand[i]);
	}
	fclose(expFile);
}
//******************************************************************************



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


  // Update SNA_buffer
  for(int s = s_l-1; s>0; s--) data->SNA_buffer[s] = data->SNA_buffer[s-1];
  data->SNA_buffer[0] = data->SNA;

  // Update PNA_buffer
  for(int p = p_l-1; p>0; p--) data->PNA_buffer[p] = data->PNA_buffer[p-1];
  data->PNA_buffer[0] = data->PNA;

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


			// EQ 2 from Scarsoglio et al. 2014
			data->Esys_lv 			= 0.59 * (data->RR[1]/data->RR[2]) + 0.91;
		}
		// HR = data->p_ursino[148+19];// + deltaHR; // eq. 18 of paper ABC.

  	data->Ts             		= 0.37	*sqrt(data->RR[0]);
  	data->Tasys          		= 0.25	*sqrt(data->RR[0]);
  	data->Tav            		= 0.19	*sqrt(data->RR[0]);
		cardiac_iter++;
  }
  // ***************************************************************************************************

  // Solver function
  CVodeSetUserData(cvode_mem, data); // you must tell cvode_mem about data. You have time dependent data.
  retval = CVode(cvode_mem, tout, y_ursino, &t, CV_NORMAL); if(check_retval(&retval, "CVode", 1)) exit(1);

// JJJ Shear Stress

// array of lengths for each flow. Might be better to choose reasonable lengths. But for now all set to Sanjay's constant parameter.
// for (i = 0; i<45; i++){
//   data->lengths[i] = data->p_ursino[101];
// }
//
// // mu is p[100].
// for(i = 1; i<45; i++){
//   data->radius[i]    	= pow(( (8.0 * data->lengths[i]  * data->p_ursino[100] )/ (M_PI * data->p_R[i]) ), 0.25 );
//   data->tau[i]       	= 4.0 * data->p_ursino[100] * (*(data->flows[i]))  / (M_PI * pow(data->radius[i],3));
// }

// *****************************************************************************
// ***** Transient outputs *****************************************************
// *****************************************************************************
  if (headersPrinted == 0) fprintf(outputInfoFile, "Output values:\n");

  fprintf(stateFile, "%f", tout);
  if (headersPrinted == 0) fprintf(outputInfoFile, "%s", "time");

  fprintf(stateFile, "\t%f", Ith(y_ursino, 10));
  if (headersPrinted == 0) fprintf(outputInfoFile, "%s", "P_a");

  fprintf(stateFile, "\t%f", data->q_ml);
  if (headersPrinted == 0) fprintf(outputInfoFile, "%s", "q_ml");

  fprintf(stateFile, "\t%f", data->q_al);
  if (headersPrinted == 0) fprintf(outputInfoFile, "%s", "q_al");

  fprintf(stateFile, "\t%f", data->q_pl);
  if (headersPrinted == 0) fprintf(outputInfoFile, "%s", "q_pl");

  fprintf(stateFile, "\t%f", data->q_mr);
  if (headersPrinted == 0) fprintf(outputInfoFile, "%s", "q_mr");

  fprintf(stateFile, "\t%f", data->q_ar);
  if (headersPrinted == 0) fprintf(outputInfoFile, "%s", "q_ar");

  fprintf(stateFile, "\t%f", data->q_pr);
  if (headersPrinted == 0) fprintf(outputInfoFile, "%s", "q_pr");

  fprintf(stateFile,"\n");
  if (headersPrinted == 0) headersPrinted = 1;
// *****************************************************************************

  tout = tout + DELTAT;
} // end of time loop.



  fclose(stateFile);
  fclose(outputInfoFile);
  fclose(postprocessedFile);

  N_VDestroy_Serial(y_ursino);
	// emxDestroyArray_real_T(X);
  CVodeFree(&cvode_mem);
  free(data);
	printf("done run %d\n",atoi(argv[1]));

	return 0; // the main must return a success exit code to the middleware.
}
