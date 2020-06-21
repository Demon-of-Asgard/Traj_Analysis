 abs_dat get_L_abs_vol(double time, char species_name[]){
	int i, j;
	int count = 0;
	double Nlumin_provided_time,avE_provided_time,tmp;
	FILE *fr;

	
	dimf dim = shape_of_file(dd2_lumins_abs_vol_fname);
	int nrows_lumins = dim.rows;
	int ncols_lumins = dim.cols;


	double times_lumins[nrows_lumins],t[nrows_lumins];
	double lumins[nrows_lumins],avE[nrows_lumins];
	
	fr = fopen(dd2_lumins_abs_vol_fname,"r");
	if(fr == NULL){
		printf("Unable to open %s, Exiting.\n",dd2_lumins_abs_vol_fname);
		exit(0);
	}
	count = 0;
	if(strcmp(species_name,"enu") == 0){
		//printf("%s\n", species_name);
		while((fscanf(fr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",\
		&t[count],&tmp,&tmp,&tmp,&lumins[count],&tmp,&tmp,&avE[count],&tmp,&tmp,&tmp,&tmp))!=EOF){
			count++;
		}
	}
	if(strcmp(species_name,"aenu") == 0){
		//printf("%s\n", species_name);
		while((fscanf(fr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",\
		&t[count],&tmp,&tmp,&tmp,&tmp,&lumins[count],&tmp,&tmp,&avE[count],&tmp,&tmp,&tmp))!=EOF){
			count++;
		}
	}
	if(strcmp(species_name,"xnu") == 0){
		//printf("%s\n", species_name);
		while((fscanf(fr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",\
		&t[count],&tmp,&tmp,&tmp,&tmp,&tmp,&lumins[count],&tmp,&tmp,&avE[count],&tmp,&tmp))!=EOF){
			lumins[count] = lumins[count]/4.0;
			count++;
		}
	}
	fclose(fr);

	/*
	 * Converting time.
	 */
	for(i=0; i<nrows_lumins-1; i++){
		times_lumins[i] = (t[i]-3441)*0.004926;
	}

	/*
	 * Finding the Number luminosity and average energy of x-neutrino at given time from the provided file.
	 */
	for(i=0; i<nrows_lumins-1;i++){
		if(times_lumins[i]<=time && times_lumins[i+1]>time){
			Nlumin_provided_time = lumins[i]+((lumins[i+1]-lumins[i])/(times_lumins[i+1]-times_lumins[i]))*(time-times_lumins[i]);
			avE_provided_time = avE[i]+((avE[i+1]-avE[i])/(times_lumins[i+1]-times_lumins[i]))*(time-times_lumins[i]);
			printf("%d. time: %lf\tNL_%s_provided: %le\t <E_%s>: %le\n",i, time,species_name,Nlumin_provided_time,\
				species_name, avE_provided_time);
			break;
		}
	}
	abs_dat abs_vol_dat;
	abs_vol_dat.NL_abs_vol = Nlumin_provided_time;
	abs_vol_dat.meanE_abs_vol = avE_provided_time;
	return(abs_vol_dat);
}
