/* Timothy J. Hunter
The composite model of Systemic circulation, detailed cerebral circulation, and baroreceptor control mechanism.
References:
The composite model is desciribed in this paper. In the pdfs/ directory, refs refer to references in this paper:
1) Ki Moo Lim, Sung Wook Choi, Byung Goo Min, and Eun Bo Shim.  Numerical Simulation of the Effect of Sodium Profile on Cardiovascular Response to Hemodialysis. Yonsei Med J 49(4):581 - 591, 2008.
2) Heldt 2002. Journal of Applied Physiology.

Units in this model are:
time: seconds.
flow: ml/s.
pressure: mmHg.
resistance: mmHg*s/ml
capacitance: ml/mmHg
elastance: mmhg/ml
*/
#include "CBF_driver.h"

int main(int argc, char* argv[])
{
    void* cvode_mem; // pointer to memory: the full state lives here.
    realtype t, tout;
    int retval;
    char* str;
    FILE *output_file, *input_data_file;

    // Create serial vector of length NEQ for I.C. and abstol
    N_Vector y_cbf = N_VNew_Serial(NEQ); // allocated memory to pointer.

    // Create userData pointer
    UserData data = (UserData)malloc(sizeof *data); // instance pointer.

    // Blood pressure states
    Ith(y_cbf, 0 + 1) = 25.7; // % Pup, 	units: mmHg
    Ith(y_cbf, 1 + 1) = 25.0; // % Pk, 		units: mmHg.=
    Ith(y_cbf, 2 + 1) = 25.0; // % Psp, 	units: mmHg
    Ith(y_cbf, 3 + 1) = 25.0; // % Pll, 	units: mmHg
    Ith(y_cbf, 4 + 1) = 5.2; // % Pab, 	units: mmHg

    // Thoracic pressure (pressure in thoracic cavity NOT blood pressure)
    Ith(y_cbf, 5 + 1) = -5.0; // % Pth, 	units: mmHg.

    // Ventricular states (BP and compliances)
    Ith(y_cbf, 6 + 1) = 10.0; // % Cl, 	units: ml/mmHg
    Ith(y_cbf, 7 + 1) = 20.0; // % Cr, 		units: ml/mmHg
    Ith(y_cbf, 8 + 1) = 12.78; // % Pl, left ventricular pressure, units: mmHg

    // More blood pressure states
    Ith(y_cbf, 9 + 1) = 80.0; // % Pa, aortic pressure, units: mmHg.
    Ith(y_cbf, 10 + 1) = 5.0; // % Psup, 		units: mmHg
    Ith(y_cbf, 11 + 1) = 5.0; // % Pinf, 		units: mmHg
    Ith(y_cbf, 12 + 1) = 5.60; // % Pr, 			units: mmHg. right ventricle pressure.
    Ith(y_cbf, 13 + 1) = 20.66; // % Ppa,			units: mmHg. Pulmonary artery pressure.
    Ith(y_cbf, 14 + 1) = 5.99; // % Ppv, 		units: mmHg. Pulmonary vein pressure.

    // Atrial states (BP and compliances)
    Ith(y_cbf, 15 + 1) = 1.0; // Cla variable, left atrial elastance.
    Ith(y_cbf, 16 + 1) = 2.0; //  Cra variable, right atrial elastance.
    Ith(y_cbf, 17 + 1) = 8.0; //  left atrial pressure initial condition, mmHg.
    Ith(y_cbf, 18 + 1) = 1.0; //  right atrial pressure initial condition, mmHg.

    // Cerebral states
    Ith(y_cbf, 19 + 1) = 9.5; // eq 1: dP_ic
    Ith(y_cbf, 20 + 1) = 25.0; // P_c
    Ith(y_cbf, 21 + 1) = 14.0; // eq. 11 : dP_v
    Ith(y_cbf, 22 + 1) = 35.5; // P_djs
    Ith(y_cbf, 23 + 1) = 35.5; // P_djs
    Ith(y_cbf, 24 + 1) = 35.5; // P_djs
    Ith(y_cbf, 25 + 1) = 35.5; // P_djs
    Ith(y_cbf, 26 + 1) = 35.5; // P_djs
    Ith(y_cbf, 27 + 1) = 35.5; // P_djs
    Ith(y_cbf, 28 + 1) = 75.5; // P_ICAl
    Ith(y_cbf, 29 + 1) = 75.5; // P_ICAr
    Ith(y_cbf, 30 + 1) = 75.5; // P_BA

    for (int i = 31; i <= 42; i++)
        Ith(y_cbf, i + 1) = 0.0; // X_aut_djs, X_co2_djs. Cerebral autoregulation states

    for (int i = 43; i <= 48; i++)
        Ith(y_cbf, i + 1) = 0.004; // C_djs. Cerebral compliance states

    // Baroreflex states
    Ith(y_cbf, 49 + 1) = 90.0; // P_aff
    Ith(y_cbf, 50 + 1) = 90.0; // X0 for P_error
    Ith(y_cbf, 51 + 1) = 0.0; // deltaHR_s
    Ith(y_cbf, 52 + 1) = 0.0; // deltaHR_v
    Ith(y_cbf, 53 + 1) = 1.0; // sigma_lv
    Ith(y_cbf, 54 + 1) = 1.0; // sigma_rv
    Ith(y_cbf, 55 + 1) = 0.0; // sigma_V
    Ith(y_cbf, 56 + 1) = 1.0; // sigma_R

    // Solver setup
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    CVodeInit(cvode_mem, RHS, 0.0, y_cbf);
    CVodeSStolerances(cvode_mem, ATOL, RTOL);
    CVDense(cvode_mem, NEQ);
    CVodeSetMaxStep(cvode_mem, DELTAT);

    // Create parameter object
    parameters(data);

    //******************************************************************************
    //*** Input/output setup *******************************************************
    //******************************************************************************

    str = malloc(32 * sizeof(char));
    sprintf(str, "statesOutput.%05d.dat", atoi(argv[1]));
    output_file = fopen(str, "w");
    free(str);

    str = malloc(32 * sizeof(char));
    sprintf(str, "input/input.%05d.dat", atoi(argv[1]));
    input_data_file = fopen(str, "r");
    free(str);

    double read_t = 0.0;

    // JJ TH Nov 3
    // calculate length of each buffer
    int s_l = (int)(data->p_ursino[37] / DELTAT);
    int p_l = (int)(data->p_ursino[26] / DELTAT);

    double deltaHR;

    tout = DELTAT;
    int cardiac_iter = 0;

    // TIME LOOP STARTS HERE
    while (tout < SIM_TIME) {

        if (tout > read_t) {
            fscanf(input_data_file, "%f %f %f", &read_t, &data->P_a, &data->PetCO2);
        }

        data->p_ursino[1] = tout; // pursino1 is time in seconds. time dependent paramter.

        // Update SNA_buffer
        for (int s = s_l - 1; s > 0; s--)
            data->SNA_buffer[s] = data->SNA_buffer[s - 1];
        data->SNA_buffer[0] = data->SNA;

        // Update PNA_buffer
        for (int p = p_l - 1; p > 0; p--)
            data->PNA_buffer[p] = data->PNA_buffer[p - 1];
        data->PNA_buffer[0] = data->PNA;

        // Cardiac interval calculations
        if (data->p_ursino[1] >= data->p_ursino[4] + data->RR[0]) {

            // fprintf(endDiastolicFile, "%f\n", data->p_ursino[1]);

            data->p_ursino[4] = data->p_ursino[1];
            data->RR[2] = data->RR[1];
            data->RR[1] = data->RR[0];

            deltaHR = (19.64 * Ith(y_cbf, 51 + 1)) - (17.95 * Ith(y_cbf, 52 + 1)) - (1.225 * pow(Ith(y_cbf, 51 + 1), 2)) + (1.357 * pow(Ith(y_cbf, 52 + 1), 2)) - (1.523 * Ith(y_cbf, 51 + 1) * Ith(y_cbf, 52 + 1));
            data->RR[0] = 60.0 / (data->p_ursino[41] + deltaHR);

            if (isAFib < 1) {
                // data->RR[0] = data->RR[0] + (data->p_ursino[43] * (data->RR[0]) / 0.033);
            } else {
                // The following is AF condition 1 and 3 as described by Scarsoglio et al. 2014
                data->p_ursino[45] = 1.0 / (-9.2 * data->RR[0] + 14.6); // r4_exponential_sample(1.0/(-9.2*data->RR[0] + 14.6));
                data->p_ursino[44] = (data->RR[0] - 1.0 / (-9.2 * data->RR[0] + 14.6)) + (data->p_ursino[43] * (data->RR[0] - 1.0 / (-9.2 * data->RR[0] + 14.6)) / 0.033);

                data->RR[0] = data->p_ursino[45] + data->p_ursino[44];

                // EQ 2 from Scarsoglio et al. 2014
                data->Esys_lv = 0.59 * (data->RR[1] / data->RR[2]) + 0.91;
            }
            // HR = data->p_ursino[148+19];// + deltaHR; // eq. 18 of paper ABC.

            data->Ts = 0.37 * sqrt(data->RR[0]);
            data->Tasys = 0.25 * sqrt(data->RR[0]);
            data->Tav = 0.19 * sqrt(data->RR[0]);
            cardiac_iter++;
        }
        // ***************************************************************************************************

        // Solver function
        CVodeSetUserData(cvode_mem, data); // you must tell cvode_mem about data. You have time dependent data.
        retval = CVode(cvode_mem, tout, y_cbf, &t, CV_NORMAL);
        if (check_retval(&retval, "CVode", 1))
            exit(1);

        // *****************************************************************************
        // ***** Transient outputs *****************************************************
        // *****************************************************************************

        fprintf(output_file, "%f", tout);

        fprintf(output_file, "\t%f", Ith(y_cbf, 10));

        // fprintf(output_file, "\t%f", data->q_ml);

        // fprintf(output_file, "\t%f", data->q_al);

        // fprintf(output_file, "\t%f", data->q_pl);

        // fprintf(output_file, "\t%f", data->q_mr);

        // fprintf(output_file, "\t%f", data->q_ar);

        // fprintf(output_file, "\t%f", data->q_pr);

        // fprintf(output_file, "\t%f", Ith(y_cbf, 50));
        // fprintf(output_file, "\t%f", Ith(y_cbf, 51));
        fprintf(output_file, "\t%f", Ith(y_cbf, 52));
        fprintf(output_file, "\t%f", Ith(y_cbf, 53));
        // fprintf(output_file, "\t%f", Ith(y_cbf, 54));
        // fprintf(output_file, "\t%f", Ith(y_cbf, 55));
        // fprintf(output_file, "\t%f", Ith(y_cbf, 56));
        // fprintf(output_file, "\t%f", Ith(y_cbf, 57));
        fprintf(output_file, "\t%f", deltaHR);
        fprintf(output_file, "\t%f", data->SNA);
        fprintf(output_file, "\t%f", data->PNA);

        fprintf(output_file, "\n");

        tout = tout + DELTAT;
    } // end of time loop.

    fclose(output_file);

    N_VDestroy_Serial(y_cbf);
    // emxDestroyArray_real_T(X);
    CVodeFree(&cvode_mem);
    free(data);
    printf("done run %d\n", atoi(argv[1]));

    return 0; // the main must return a success exit code to the middleware.
}
