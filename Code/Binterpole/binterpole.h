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

	double Nnu50Theta1_t1, Nnu50Theta1_t2, Nnu50Theta1_t, Enu50Theta1_t;
	double Nnu50Theta2_t1, Nnu50Theta2_t2, Nnu50Theta2_t, Enu50Theta2_t;
	double Nnu50Theta_t,Enu50Theta_t;

	double Nnu100Theta1_t1, Nnu100Theta1_t2, Nnu100Theta1_t, Enu100Theta1_t;
	double Nnu100Theta2_t1, Nnu100Theta2_t2, Nnu100Theta2_t, Enu100Theta2_t;
	double Nnu100Theta_t, Enu100Theta_t;

	double NnuTheta_rt, EnuTheta_rt;



	double Nanu50Theta1_t1, Nanu50Theta1_t2, Nanu50Theta1_t, Eanu50Theta1_t;
	double Nanu50Theta2_t1, Nanu50Theta2_t2, Nanu50Theta2_t, Eanu50Theta2_t;
	double Nanu50Theta_t, Eanu50Theta_t;

	double Nanu100Theta1_t1, Nanu100Theta1_t2, Nanu100Theta1_t, Eanu100Theta1_t;
	double Nanu100Theta2_t1, Nanu100Theta2_t2, Nanu100Theta2_t, Eanu100Theta2_t;
	double Nanu100Theta_t, Eanu100Theta_t;

	double NanuTheta_rt, EanuTheta_rt;




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
	theta = acos(z/r)*180.0/M_PI;// To deg.

	int bin_number 	 = 0;  //1 + (int) (theta*(180.0/M_PI))/10;
	double thetabin1 = 0.0;//((bin_number-1)*10.0)*M_PI/180.0;
	double thetabin2 = 0.0;//((bin_number)*10.0)*M_PI/180.0;

	if(theta>=0.0 && theta<=10.0){bin_number = 1; thetabin1  =  0.0*M_PI/180.0; thetabin2  = 10.0*M_PI/180.0;}
	if(theta>10.0 && theta<=20.0){bin_number = 2; thetabin1  = 10.0*M_PI/180.0; thetabin2  = 20.0*M_PI/180.0;}
	if(theta>20.0 && theta<=30.0){bin_number = 3; thetabin1  = 20.0*M_PI/180.0; thetabin2  = 30.0*M_PI/180.0;}
	if(theta>30.0 && theta<=40.0){bin_number = 4; thetabin1  = 30.0*M_PI/180.0; thetabin2  = 40.0*M_PI/180.0;}
	if(theta>40.0 && theta<=50.0){bin_number = 5; thetabin1  = 40.0*M_PI/180.0; thetabin2  = 50.0*M_PI/180.0;}
	if(theta>50.0 && theta<=60.0){bin_number = 6; thetabin1  = 50.0*M_PI/180.0; thetabin2  = 60.0*M_PI/180.0;}
	if(theta>60.0 && theta<=70.0){bin_number = 7; thetabin1  = 60.0*M_PI/180.0; thetabin2  = 70.0*M_PI/180.0;}
	if(theta>70.0 && theta<=80.0){bin_number = 8; thetabin1  = 70.0*M_PI/180.0; thetabin2  = 80.0*M_PI/180.0;}
	if(theta>80.0 && theta<=90.0){bin_number = 9; thetabin1  = 80.0*M_PI/180.0; thetabin2  = 90.0*M_PI/180.0;}
	if(theta>90.0 && theta<=100.0){bin_number = 10; thetabin1  = 90.0*M_PI/180.0; thetabin2  = 100.0*M_PI/180.0;}
	if(theta>100.0 && theta<=110.0){bin_number = 11; thetabin1  = 100.0*M_PI/180.0; thetabin2  = 110.0*M_PI/180.0;}
	if(theta>110.0 && theta<=120.0){bin_number = 12; thetabin1  = 110.0*M_PI/180.0; thetabin2  = 120.0*M_PI/180.0;}
	if(theta>120.0 && theta<=130.0){bin_number = 13; thetabin1  = 120.0*M_PI/180.0; thetabin2  = 130.0*M_PI/180.0;}
	if(theta>130.0 && theta<=140.0){bin_number = 14; thetabin1  = 130.0*M_PI/180.0; thetabin2  = 140.0*M_PI/180.0;}
	if(theta>140.0 && theta<=150.0){bin_number = 15; thetabin1  = 140.0*M_PI/180.0; thetabin2  = 150.0*M_PI/180.0;}
	if(theta>150.0 && theta<=160.0){bin_number = 16; thetabin1  = 150.0*M_PI/180.0; thetabin2  = 160.0*M_PI/180.0;}
	if(theta>160.0 && theta<=170.0){bin_number = 17; thetabin1  = 160.0*M_PI/180.0; thetabin2  = 170.0*M_PI/180.0;}
	if(theta>170.0 && theta<=180.0){bin_number = 18; thetabin1  = 170.0*M_PI/180.0; thetabin2  = 180.0*M_PI/180.0;}

	theta = theta*M_PI/180.0; // Back to rad.

	int row_number = 0;
	for(i=0;Lnu50[i][0]<time;i++){
		row_number = i;
	}

	double t1 = Lnu50[row_number][0];
	double t2 = Lnu50[row_number+1][0];

	printf("theta: %le  Bin#: %d\n", theta, bin_number);
	printf("row_number = %d,t1=%le\t t2=%le\n\v", row_number,t1,t2);

	double r50 = 50.0e5;
	double r100 = 100.0e5;
	double Delta_A50 = 2.0*M_PI*r50*r50*fabs(cos(thetabin1)-cos(thetabin2));
	double Delta_A100 = 2.0*M_PI*r100*r100*fabs(cos(thetabin1)-cos(thetabin2));

	//==================================================\nu==========================================================================
	//===============================================================================================================================
	if(bin_number==18){
		//===============================================================================================================================

		Nnu50Theta1_t1 = Lnu50[row_number][bin_number]/(c*Delta_A50*avEnu50[row_number][bin_number]);
		Nnu50Theta1_t2 = Lnu50[row_number+1][bin_number]/(c*Delta_A50*avEnu50[row_number+1][bin_number]);
		Nnu50Theta1_t  = NewtonsInterpolation(time, t1, t2, Nnu50Theta1_t1, Nnu50Theta1_t2);
		Enu50Theta1_t  = NewtonsInterpolation(time, t1, t2, avEnu50[row_number][bin_number], avEnu50[row_number+1][bin_number]);

		Nnu50Theta_t = Nnu50Theta1_t;
		Enu50Theta_t = Enu50Theta1_t;
		printf("Nnu50Theta_t: %le\n",Nnu50Theta_t);
		printf("Enu50Theta_t: %le\n",Enu50Theta_t);

		//==================================================================================================================================

		Nnu100Theta1_t1 = Lnu100[row_number][bin_number]/(c*Delta_A100*avEnu100[row_number][bin_number]);
		Nnu100Theta1_t2 = Lnu100[row_number+1][bin_number]/(c*Delta_A100*avEnu100[row_number+1][bin_number]);
		Nnu100Theta1_t  = NewtonsInterpolation(time, t1, t2, Nnu100Theta1_t1, Nnu100Theta1_t2);
		Enu100Theta1_t  = NewtonsInterpolation(time, t1, t2, avEnu100[row_number][bin_number], avEnu100[row_number+1][bin_number]);

		Nnu100Theta_t = Nnu100Theta1_t;
		Enu100Theta_t = Enu100Theta1_t;

		printf("Nnu100Theta_t: %le\n",Nnu100Theta_t);
		printf("Enu100Theta_t: %le\n",Enu100Theta_t);

		//===================================================================================================================================
		NnuTheta_rt = LogLinInterpolation(r,r50,r100,Nnu50Theta_t, Nnu100Theta_t);
		EnuTheta_rt = NewtonsInterpolation(r, r50, r100, Enu50Theta_t, Enu100Theta_t);

		printf("NnuTheta_rt: %le\n",NnuTheta_rt);
		printf("EnuTheta_rt: %le\n",EnuTheta_rt);
		//====================================================\bar\nu=======================================================================
		//==================================================================================================================================

		Nanu50Theta1_t1 = Lanu50[row_number][bin_number]/(c*Delta_A50*avEanu50[row_number][bin_number]);
		Nanu50Theta1_t2 = Lanu50[row_number+1][bin_number]/(c*Delta_A50*avEanu50[row_number+1][bin_number]);
		Nanu50Theta1_t  = NewtonsInterpolation(time, t1, t2, Nanu50Theta1_t1, Nanu50Theta1_t2);
		Eanu50Theta1_t  = NewtonsInterpolation(time, t1, t2, avEanu50[row_number][bin_number], avEanu50[row_number+1][bin_number]);

		Nanu50Theta_t = Nanu50Theta1_t;
		Eanu50Theta_t = Eanu50Theta1_t;
		//=================================================================================================================================

		Nanu100Theta1_t1 = Lanu100[row_number][bin_number]/(c*Delta_A100*avEanu100[row_number][bin_number]);
		Nanu100Theta1_t2 = Lanu100[row_number+1][bin_number]/(c*Delta_A100*avEanu100[row_number+1][bin_number]);
		Nanu100Theta1_t  = NewtonsInterpolation(time, t1, t2, Nanu100Theta1_t1, Nanu100Theta1_t2);
		Eanu100Theta1_t  = NewtonsInterpolation(time, t1, t2, avEanu100[row_number][bin_number], avEanu100[row_number+1][bin_number]);

		Nanu100Theta_t = Nanu100Theta1_t;
		Eanu100Theta_t = Eanu100Theta1_t;
		//==================================================================================================================================

		NanuTheta_rt = LogLinInterpolation(r,r50,r100,Nanu50Theta_t, Nanu100Theta_t);
		EanuTheta_rt = NewtonsInterpolation(r, r50, r100, Eanu50Theta_t, Eanu100Theta_t);
		//===================================================================================================================================
	}
	if(bin_number<18){
		//===============================================================================================================================

		Nnu50Theta1_t1 = Lnu50[row_number][bin_number]/(c*Delta_A50*avEnu50[row_number][bin_number]);
		Nnu50Theta1_t2 = Lnu50[row_number+1][bin_number]/(c*Delta_A50*avEnu50[row_number+1][bin_number]);
		Nnu50Theta1_t  = NewtonsInterpolation(time, t1, t2, Nnu50Theta1_t1, Nnu50Theta1_t2);
		Enu50Theta1_t  = NewtonsInterpolation(time, t1, t2, avEnu50[row_number][bin_number], avEnu50[row_number+1][bin_number]);

		Nnu50Theta2_t1 = Lnu50[row_number][bin_number+1]/(c*Delta_A50*avEnu50[row_number][bin_number+1]);
		Nnu50Theta2_t2 = Lnu50[row_number+1][bin_number+1]/(c*Delta_A50*avEnu50[row_number+1][bin_number+1]);
		Nnu50Theta2_t  = NewtonsInterpolation(time, t1, t2, Nnu50Theta2_t1, Nnu50Theta2_t2);
		Enu50Theta2_t  = NewtonsInterpolation(time, t1, t2, avEnu50[row_number][bin_number+1], avEnu50[row_number+1][bin_number+1]);

		Nnu50Theta_t = NewtonsInterpolation(theta, thetabin1, thetabin2, Nnu50Theta1_t, Nnu50Theta2_t);
		Enu50Theta_t = NewtonsInterpolation(theta, thetabin1, thetabin2, Enu50Theta1_t, Enu50Theta2_t);


		//==================================================================================================================================

		Nnu100Theta1_t1 = Lnu100[row_number][bin_number]/(c*Delta_A100*avEnu100[row_number][bin_number]);
		Nnu100Theta1_t2 = Lnu100[row_number+1][bin_number]/(c*Delta_A100*avEnu100[row_number+1][bin_number]);
		Nnu100Theta1_t  = NewtonsInterpolation(time, t1, t2, Nnu100Theta1_t1, Nnu100Theta1_t2);
		Enu100Theta1_t  = NewtonsInterpolation(time, t1, t2, avEnu100[row_number][bin_number], avEnu100[row_number+1][bin_number]);

		Nnu100Theta2_t1 = Lnu100[row_number][bin_number+1]/(c*Delta_A100*avEnu100[row_number][bin_number+1]);
		Nnu100Theta2_t2 = Lnu100[row_number+1][bin_number+1]/(c*Delta_A100*avEnu100[row_number+1][bin_number+1]);
		Nnu100Theta2_t  = NewtonsInterpolation(time, t1, t2, Nnu100Theta2_t1, Nnu100Theta2_t2);
		Enu100Theta2_t  = NewtonsInterpolation(time, t1, t2, avEnu100[row_number][bin_number+1], avEnu100[row_number+1][bin_number+1]);

		Nnu100Theta_t = NewtonsInterpolation(theta, thetabin1, thetabin2, Nnu100Theta1_t, Nnu100Theta2_t);
		Enu100Theta_t = NewtonsInterpolation(theta, thetabin1, thetabin2, Enu100Theta1_t, Enu100Theta2_t);

		//===================================================================================================================================

		NnuTheta_rt = LogLinInterpolation(r,r50,r100,Nnu50Theta_t, Nnu100Theta_t);
		EnuTheta_rt = NewtonsInterpolation(r, r50, r100, Enu50Theta_t, Enu100Theta_t);
		//double NnuTheta_rt = LinInverseR2Interpolation(r,r50,r100,Nnu50Theta_t, Nnu100Theta_t);

		//====================================================\bar\nu=======================================================================
		//==================================================================================================================================

		Nanu50Theta1_t1 = Lanu50[row_number][bin_number]/(c*Delta_A50*avEanu50[row_number][bin_number]);
		Nanu50Theta1_t2 = Lanu50[row_number+1][bin_number]/(c*Delta_A50*avEanu50[row_number+1][bin_number]);
		Nanu50Theta1_t  = NewtonsInterpolation(time, t1, t2, Nanu50Theta1_t1, Nanu50Theta1_t2);
		Eanu50Theta1_t  = NewtonsInterpolation(time, t1, t2, avEanu50[row_number][bin_number], avEanu50[row_number+1][bin_number]);

		Nanu50Theta2_t1 = Lanu50[row_number][bin_number+1]/(c*Delta_A50*avEanu50[row_number][bin_number+1]);
		Nanu50Theta2_t2 = Lanu50[row_number+1][bin_number+1]/(c*Delta_A50*avEanu50[row_number+1][bin_number+1]);
		Nanu50Theta2_t  = NewtonsInterpolation(time, t1, t2, Nanu50Theta2_t1, Nanu50Theta2_t2);
		Eanu50Theta2_t  = NewtonsInterpolation(time, t1, t2, avEanu50[row_number][bin_number+1], avEanu50[row_number+1][bin_number+1]);

		Nanu50Theta_t = NewtonsInterpolation(theta, thetabin1, thetabin2, Nanu50Theta1_t, Nanu50Theta2_t);
		Eanu50Theta_t = NewtonsInterpolation(theta, thetabin1, thetabin2, Eanu50Theta1_t, Eanu50Theta2_t);

		//=================================================================================================================================

		Nanu100Theta1_t1 = Lanu100[row_number][bin_number]/(c*Delta_A100*avEanu100[row_number][bin_number]);
		Nanu100Theta1_t2 = Lanu100[row_number+1][bin_number]/(c*Delta_A100*avEanu100[row_number+1][bin_number]);
		Nanu100Theta1_t  = NewtonsInterpolation(time, t1, t2, Nanu100Theta1_t1, Nanu100Theta1_t2);
		Eanu100Theta1_t  = NewtonsInterpolation(time, t1, t2, avEanu100[row_number][bin_number], avEanu100[row_number+1][bin_number]);

		Nanu100Theta2_t1 = Lanu100[row_number][bin_number+1]/(c*Delta_A100*avEanu100[row_number][bin_number+1]);
		Nanu100Theta2_t2 = Lanu100[row_number+1][bin_number+1]/(c*Delta_A100*avEanu100[row_number+1][bin_number+1]);
		Nanu100Theta2_t  = NewtonsInterpolation(time, t1, t2, Nanu100Theta2_t1, Nanu100Theta2_t2);
		Eanu100Theta2_t  = NewtonsInterpolation(time, t1, t2, avEanu100[row_number][bin_number+1], avEanu100[row_number+1][bin_number]);

		Nanu100Theta_t = NewtonsInterpolation(theta, thetabin1, thetabin2, Nanu100Theta1_t, Nanu100Theta2_t);
		Eanu100Theta_t = NewtonsInterpolation(theta, thetabin1, thetabin2, Eanu100Theta1_t, Eanu100Theta2_t);

		//==================================================================================================================================

		NanuTheta_rt = LogLinInterpolation(r,r50,r100,Nanu50Theta_t, Nanu100Theta_t);
		EanuTheta_rt = NewtonsInterpolation(r, r50, r100, Eanu50Theta_t, Eanu100Theta_t);
		//double NanuTheta_rt = LinInverseR2Interpolation(r,r50,r100,Nanu50Theta_t, Nanu100Theta_t);

		//===================================================================================================================================

	}
	bin_polated binout;

	binout.Nnu_rt    = NnuTheta_rt;
	binout.avEnu_rt  = EnuTheta_rt/1.60218e-6;
	binout.Nanu_rt   = NanuTheta_rt;
	binout.avEanu_rt = EanuTheta_rt/1.60218e-6;
	printf("theta: %le bin_number: %d\n", theta*180.0/M_PI,bin_number);
	printf("Nnu_rt: %le Nanu_rt: %le\n",NnuTheta_rt,NanuTheta_rt);
	printf("Enu_rt: %le Eanu_rt: %le\n",EnuTheta_rt/1.60218e-6,EanuTheta_rt/1.60218e-6);

	//===================================================================================================================================
	printf("\v");

	return(binout);
}
