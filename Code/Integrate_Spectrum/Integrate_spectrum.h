values_at_loc Number_density_and_avEnergy(char *prof_fname){
	dimf dim = shape_of_file(prof_fname);
	int nrows = dim.rows;
	int ncols = dim.cols;

	FILE *fr, *fw;
	int i, j, count = 0, N = 0;
	double du = ((double)(2.0/NU));
	double dp = ((double)(2.0*M_PI/NP));
	double tmp, nu_data[nrows], Enu_data[nrows];
	double nu_prof[NU][NP], Enu_prof[NU][NP];
	double Nnu = 0.0,avEnu = 0.0;
	values_at_loc values;

	fr = fopen(prof_fname, "r");
	count = 0;
	while((fscanf(fr,"%lf%lf%lf%lf",&tmp,&tmp,&nu_data[count],&Enu_data[count]))!=EOF){
		count++;
	}
	fclose(fr);
	for(i=0; i<NU; i++){
		for(j=0; j<NP; j++){
			N = i*NP+j;
			nu_prof[i][j] = nu_data[N];
			Enu_prof[i][j] = Enu_data[N];
		}
	}

	for(i=0; i<NU; i++){
		for(j=0; j<NP; j++){
			Nnu += (nu_prof[i][j]*du*dp);
			avEnu += (Enu_prof[i][j]*du*dp);
		}
	}
	avEnu = avEnu/Nnu;

	values.n_nu = Nnu;
	values.av_e_nu = avEnu;
	//printf("#xnu @ given loc: %le, <E_xnu> @ given loc: %le\n", *values_at_loc, *(values_at_loc + 1));
	return(values);
}
