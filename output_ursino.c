// *****************************************************************************************************************
// my OUTPUTS


// if (cardiac_iter >= 20){
	// START OF CARDIAC OUTPUT
	//
	// if((data->p_ursino[105] >= DELTAT)&& (data->p_ursino[105]<data->p_ursino[11])) {
	// 		cardiac_output = cardiac_output + (data->Qlo) * DELTAT;
	//     flowpermin[0] += data->CBF*DELTAT;
	//     flowpermin[1] += data->q_ml*DELTAT;
	//     flowpermin[2] += data->q_al*DELTAT;
	//     flowpermin[3] += data->q_pl*DELTAT;
	// 		flowpermin[4] += data->q_ml*DELTAT;
	//     flowpermin[5] += data->q_al*DELTAT;
	//     flowpermin[6] += data->q_pl*DELTAT;
	//     // printf("%f\t%f\t%f\t%f\n", data->p_ursino[105], data->p_ursino[11],  cardiac_output);
	// } else if((data->p_ursino[105] < DELTAT) && (cardiac_output > 0.0)) {
	// 	// printf("RR is %f",data->RR[0]);
	// 		fprintf(cardiac_output_file, "%f\t%f\t%f\t%f\t%f\t", tout, data->RR[0], data->p_ursino[170], data->p_ursino[171], cardiac_output*(60.0/data->RR[0]));
	//     for(int f = 0; f<7; f++){
	//       fprintf(cardiac_output_file, "%f\t", flowpermin[f]*(60.0/data->RR[0]));
	//       flowpermin[f] = 0.0;
	//     }
	//     fprintf(cardiac_output_file, "\n");
	// 		cardiac_output = 0.0;
	// }

	// END OF CARDIAC OUTPUT

	// can choose to only output after iout reaches a certain percentage of NOUT i.e 90%
	// if(tout>=SIMTIME-10.0){
	// fprintf(cardiac_output_file, "%f\t", tout);
	// fprintf(cardiac_output_file, "%f\t",Ith(y_ursino, 13)); // arterial pressure
	// fprintf(cardiac_output_file,"\n");



	/*
	CEREBRAL***.dat
	1) Time
	2-7) q_js
	8) Pa
	9) P_sup
	10) P_ic
	11) P_c
	12) P_v
	13-18) P_djs
	19) P_ICAl
	20) P_ICAr
	21) P_B
	*/
// if ((int)(tout*1000) % 10 == 0){
if (headersPrinted == 0){
  fprintf(cerebral, "%s\t", "time");
  fprintf(cerebral, "%s\t", "P_a");
  fprintf(cerebral, "%s\t%s\t%s\t%s\t%s\t%s\n", "q_ml", "q_al", "q_pl", "q_mr", "q_ar", "q_pr");
  fprintf(cerebral, "\n");
  headersPrinted = 1;
}
	fprintf(cerebral, "%f\t", tout);
  fprintf(cerebral, "%s\n", Ith(y_ursino, 10));

	fprintf(cerebral, "%f\t%f\t%f\t%f\t%f\t%f\t",data->q_ml, data->q_al, data->q_pl, data->q_mr, data->q_ar, data->q_pr);
	// fprintf(cerebral, "%f\t",data->q_ICAl);
	// fprintf(cerebral, "%f\t",data->temp); // currently P_PCAl
	// fprintf(cerebral, "%f\t",Ith(y_ursino, 55));
	// fprintf(cerebral, "%f\t",Ith(y_ursino, 78));
	// fprintf(cerebral, "%f\t",Ith(y_ursino, 79));
	// fprintf(cerebral, "%f\t",Ith(y_ursino, 80));
	// fprintf(cerebral, "%f\t",data->q_ICAr);
	// fprintf(cerebral, "%f\t",data->q_MCAl);
	// fprintf(cerebral, "%f\t",data->q_MCAr);
	// fprintf(cerebral, "%f\t",data->q_ACAl);
	// fprintf(cerebral, "%f\t",data->q_ACAr);
	// fprintf(cerebral, "%f\t",data->q_PCAl);
	// fprintf(cerebral, "%f\t",data->q_PCAr);
	// fprintf(cerebral, "%f\t",data->q_PCoAl);
	// fprintf(cerebral, "%f\t",data->q_PCoAr);

	fprintf(cerebral,"\n");
// }
// }
