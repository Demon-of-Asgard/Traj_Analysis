/*
 * Program to extrapolate the bin data on to the trajectory.
 */
bin_polated binterpole(double time, double x, double z){
	x = x*1.0e5; 
	z = z*1.0e5;

	int i = 0;
	int j = 0;
	double  r, theta;
	double tmp;

	FILE *frEanu50, *frEanu100, *frEnu50, *frEnu100;
	FILE *frLanu50, *frLanu100, *frLnu50, *frLnu100;
	FILE *fw;

	dimf dim = shape_of_file("../BinData/luminnu50.dat");

	double Lnu50[dim.rows][dim.cols], Lnu100[dim.rows][dim.cols], avEnu50[dim.rows][dim.cols],\
	avEnu100[dim.rows][dim.cols];

	double Lanu50[dim.rows][dim.cols], Lanu100[dim.rows][dim.cols], avEanu50[dim.rows][dim.cols],\
	avEanu100[dim.rows][dim.cols];

	frLnu50   = fopen("../BinData/luminnu50.dat", "r");
	frLnu100  = fopen("../BinData/luminnu100.dat","r");
	frLanu50  = fopen("../BinData/luminanu50.dat", "r");
	frLanu100 = fopen("../BinData/luminanu100.dat","r");

	frEnu50   = fopen("../BinData/avEnu50.dat", "r");
	frEanu50  = fopen("../BinData/avEanu50.dat", "r");
	frEnu100  = fopen("../BinData/avEnu100.dat","r");
	frEanu100 = fopen("../BinData/avEanu100.dat","r");

	for(i=0;i<dim.rows; i++){
		for(j=0; (j<dim.cols);j++){
			fscanf(frLnu50,"%lf", &Lnu50[i][j]);
			fscanf(frLnu100, "%lf", &Lnu100[i][j]);
			fscanf(frLanu50, "%lf", &Lanu50[i][j]);
			fscanf(frLanu100, "%lf", &Lanu100[i][j]);

			fscanf(frEnu50, "%lf", &avEnu50[i][j]);
			fscanf(frEanu50, "%lf", &avEanu50[i][j]);
			fscanf(frEnu100, "%lf", &avEnu100[i][j]);
			fscanf(frEanu100, "%lf", &avEanu100[i][j]);

			avEnu50[i][j]   = fabs(avEnu50[i][j])*1.60218e-6;
			avEanu50[i][j]  = fabs(avEanu50[i][j])*1.60218e-6;
			avEnu100[i][j]  = fabs(avEnu100[i][j])*1.60218e-6;
			avEanu100[i][j] = fabs(avEanu100[i][j])*1.60218e-6;
		}
	}

	fclose(frLnu50);
	fclose(frLnu100);
	fclose(frEnu50);
	fclose(frEnu100);

	fclose(frLanu50);
	fclose(frLanu100);
	fclose(frEanu50);
	fclose(frEanu100);

	/*
	 * subracting offset from time.
	 */
	double time_offset = 3441*0.004926;
	for(i=0;i<dim.rows;i++){
		Lnu50[i][0]    = Lnu50[i][0] - time_offset;
		Lnu100[i][0]   = Lnu100[i][0] - time_offset;
		avEnu50[i][0]  = avEnu50[i][0] - time_offset;
		avEnu100[i][0] = avEnu100[i][0] - time_offset;

		Lanu50[i][0]    = Lanu50[i][0] - time_offset;
		Lanu100[i][0]   = Lanu100[i][0] - time_offset;
		avEanu50[i][0]  = avEanu50[i][0] - time_offset;
		avEanu100[i][0] = avEanu100[i][0] - time_offset;
	}

	r = sqrt(x*x + z*z);
	theta = round((180.0/M_PI)*acos(z/r)) * (M_PI/180.0);
	//theta = acos(z/r);


	//printf("theta = %le\n\v",theta*(180.0/M_PI));

	int bin_number = 1 + (int) (theta*(180.0/M_PI))/10;
	double thetabin1 = ((bin_number-1)*10.0)*M_PI/180.0;
	double thetabin2 = ((bin_number)*10.0)*M_PI/180.0;


	int row_number = 0;
	for(i=0;Lnu50[i][0]<time;i++){
		row_number = i;
	}
	double t1 = Lnu50[row_number][0];
	double t2 = Lnu50[row_number+1][0];

	printf("row_number = %d,t1=%le\tt2=%le\n\v", row_number,t1,t2);

	double r50 = 50.0e5;
	double r100 = 100.0e5;
	double Delta_A50 = 2.0*M_PI*r50*r50*fabs(cos(thetabin1)-cos(thetabin2));
	double Delta_A100 = 2.0*M_PI*r100*r100*fabs(cos(thetabin1)-cos(thetabin2));

	//==================================================\nu==========================================================================
	//===============================================================================================================================

	double Nnu50Theta1_t1 = Lnu50[row_number][bin_number]/(c*Delta_A50*avEnu50[row_number][bin_number]);
	double Nnu50Theta1_t2 = Lnu50[row_number+1][bin_number]/(c*Delta_A50*avEnu50[row_number+1][bin_number]);
	double Nnu50Theta1_t  = NewtonsInterpolation(time, t1, t2, Nnu50Theta1_t1, Nnu50Theta1_t2);
	double Enu50Theta1_t  = NewtonsInterpolation(time, t1, t2, avEnu50[row_number][bin_number], avEnu50[row_number+1][bin_number]);

	double Nnu50Theta2_t1 = Lnu50[row_number][bin_number+1]/(c*Delta_A50*avEnu50[row_number][bin_number+1]);
	double Nnu50Theta2_t2 = Lnu50[row_number+1][bin_number+1]/(c*Delta_A50*avEnu50[row_number+1][bin_number+1]);
	double Nnu50Theta2_t  = NewtonsInterpolation(time, t1, t2, Nnu50Theta2_t1, Nnu50Theta2_t2);
	double Enu50Theta2_t  = NewtonsInterpolation(time, t1, t2, avEnu50[row_number][bin_number+1], avEnu50[row_number+1][bin_number+1]);

	double Nnu50Theta_t = LogLinInterpolation(theta, thetabin1, thetabin2, Nnu50Theta1_t, Nnu50Theta2_t);
	double Enu50Theta_t = LogLinInterpolation(theta, thetabin1, thetabin2, Enu50Theta1_t, Enu50Theta2_t);


	//==================================================================================================================================

	double Nnu100Theta1_t1 = Lnu100[row_number][bin_number]/(c*Delta_A100*avEnu100[row_number][bin_number]);
	double Nnu100Theta1_t2 = Lnu100[row_number+1][bin_number]/(c*Delta_A100*avEnu100[row_number+1][bin_number]);
	double Nnu100Theta1_t  = NewtonsInterpolation(time, t1, t2, Nnu100Theta1_t1, Nnu100Theta1_t2);
	double Enu100Theta1_t  = NewtonsInterpolation(time, t1, t2, avEnu100[row_number][bin_number], avEnu100[row_number+1][bin_number]);
 
	double Nnu100Theta2_t1 = Lnu100[row_number][bin_number+1]/(c*Delta_A100*avEnu100[row_number][bin_number+1]);
	double Nnu100Theta2_t2 = Lnu100[row_number+1][bin_number+1]/(c*Delta_A100*avEnu100[row_number+1][bin_number+1]);
	double Nnu100Theta2_t  = NewtonsInterpolation(time, t1, t2, Nnu100Theta2_t1, Nnu100Theta2_t2);
	double Enu100Theta2_t  = NewtonsInterpolation(time, t1, t2, avEnu100[row_number][bin_number+1], avEnu100[row_number+1][bin_number+1]);

	double Nnu100Theta_t = LogLinInterpolation(theta, thetabin1, thetabin2, Nnu100Theta1_t, Nnu100Theta2_t);
	double Enu100Theta_t = LogLinInterpolation(theta, thetabin1, thetabin2, Enu100Theta1_t, Enu100Theta2_t);

	//===================================================================================================================================

	double NnuTheta_rt = LogLinInterpolation(r,r50,r100,Nnu50Theta_t, Nnu100Theta_t);
	double EnuTheta_rt = LogLinInterpolation(r, r50, r100, Enu50Theta_t, Enu100Theta_t);
	//double NnuTheta_rt = LinInverseR2Interpolation(r,r50,r100,Nnu50Theta_t, Nnu100Theta_t);
	
	//====================================================\bar\nu=======================================================================
	//==================================================================================================================================

	double Nanu50Theta1_t1 = Lanu50[row_number][bin_number]/(c*Delta_A50*avEanu50[row_number][bin_number]);
	double Nanu50Theta1_t2 = Lanu50[row_number+1][bin_number]/(c*Delta_A50*avEanu50[row_number+1][bin_number]);
	double Nanu50Theta1_t  = NewtonsInterpolation(time, t1, t2, Nanu50Theta1_t1, Nanu50Theta1_t2);
	double Eanu50Theta1_t  = NewtonsInterpolation(time, t1, t2, avEanu50[row_number][bin_number], avEanu50[row_number+1][bin_number]);

	double Nanu50Theta2_t1 = Lanu50[row_number][bin_number+1]/(c*Delta_A50*avEanu50[row_number][bin_number+1]);
	double Nanu50Theta2_t2 = Lanu50[row_number+1][bin_number+1]/(c*Delta_A50*avEanu50[row_number+1][bin_number+1]);
	double Nanu50Theta2_t  = NewtonsInterpolation(time, t1, t2, Nanu50Theta2_t1, Nanu50Theta2_t2);
	double Eanu50Theta2_t  = NewtonsInterpolation(time, t1, t2, avEanu50[row_number][bin_number+1], avEanu50[row_number+1][bin_number+1]);

	double Nanu50Theta_t = LogLinInterpolation(theta, thetabin1, thetabin2, Nanu50Theta1_t, Nanu50Theta2_t);
	double Eanu50Theta_t = LogLinInterpolation(theta, thetabin1, thetabin2, Eanu50Theta1_t, Eanu50Theta2_t);

	//=================================================================================================================================

	double Nanu100Theta1_t1 = Lanu100[row_number][bin_number]/(c*Delta_A100*avEanu100[row_number][bin_number]);
	double Nanu100Theta1_t2 = Lanu100[row_number+1][bin_number]/(c*Delta_A100*avEanu100[row_number+1][bin_number]);
	double Nanu100Theta1_t  = NewtonsInterpolation(time, t1, t2, Nanu100Theta1_t1, Nanu100Theta1_t2);
	double Eanu100Theta1_t  = NewtonsInterpolation(time, t1, t2, avEanu100[row_number][bin_number], avEanu100[row_number+1][bin_number]);

	double Nanu100Theta2_t1 = Lanu100[row_number][bin_number+1]/(c*Delta_A100*avEanu100[row_number][bin_number+1]);
	double Nanu100Theta2_t2 = Lanu100[row_number+1][bin_number+1]/(c*Delta_A100*avEanu100[row_number+1][bin_number+1]);
	double Nanu100Theta2_t  = NewtonsInterpolation(time, t1, t2, Nanu100Theta2_t1, Nanu100Theta2_t2);
	double Eanu100Theta2_t  = NewtonsInterpolation(time, t1, t2, avEanu100[row_number][bin_number+1], avEanu100[row_number+1][bin_number]); 

	double Nanu100Theta_t = LogLinInterpolation(theta, thetabin1, thetabin2, Nanu100Theta1_t, Nanu100Theta2_t);
	double Eanu100Theta_t = LogLinInterpolation(theta, thetabin1, thetabin2, Eanu100Theta1_t, Eanu100Theta2_t);

	//==================================================================================================================================

	double NanuTheta_rt = LogLinInterpolation(r,r50,r100,Nanu50Theta_t, Nanu100Theta_t);
	double EanuTheta_rt = LogLinInterpolation(r, r50, r100, Eanu50Theta_t, Eanu100Theta_t);
	//double NanuTheta_rt = LinInverseR2Interpolation(r,r50,r100,Nanu50Theta_t, Nanu100Theta_t);
	
	//===================================================================================================================================

	bin_polated binout;

	binout.Nnu_rt    = NnuTheta_rt;
	binout.avEnu_rt  = EnuTheta_rt/1.60218e-6;
	binout.Nanu_rt   = NanuTheta_rt;
	binout.avEanu_rt = EanuTheta_rt/1.60218e-6;

	//===================================================================================================================================
	printf("\v");

	return(binout);
}