char* NormalizeNFlux(char species_name[], char *density_on_surface_fname, char *temperature_on_surface_fname, char *mu_on_surface_fname, double time){
	//printf("Inside NormFlux. \n");

	FILE *fr, *fw, *fw2;
	char *outfname_avE_surf;
	dimf dim;
	int nrows_surf, ncols_surf,nrows_lumins,ncols_lumins, i, j, k,N,count;
	double x0,z0,x1,z1,dl,r0,r1,av_r, av_z, theta,sintheta;
	double av_dens,dA;
	double Nlumin_estimated_time,Nlumin_provided_time,avE_provided_time, tmp, Nlumins_time;
	double flux_renormalization_factor;
	double mubar;

	/*
	Finding the #rows and #columns in the file 'density_on_surface_fname'.
	*/
	dim = shape_of_file(density_on_surface_fname);
	nrows_surf = dim.rows;
	ncols_surf = dim.cols;
	//printf("--------------->%d\t%d\t in %s\n",*dim,*(dim+1),density_on_surface_fname);

	double x[nrows_surf],z[nrows_surf];
	double xT[nrows_surf],zT[nrows_surf];
	double flux_xz[nrows_surf];
	double density[nrows_surf];
	double T_on_surface[nrows_surf];
	double mu_on_surface[nrows_surf];
	double avE_surface[nrows_surf];

	fr = fopen(density_on_surface_fname,"r");
	if(fr == NULL){
		printf("Unable to open %s. exiting.\n", density_on_surface_fname);
		exit(0);
	}
	else{
		count = 0;
		while((fscanf(fr,"%lf%lf%lf",&x[count], &z[count], &density[count]))!=EOF){
			count++;
		}
	}
	fclose(fr);

	fr = fopen(temperature_on_surface_fname,"r");
	if(fr == NULL){
		printf("Unable to open %s. exiting.\n", temperature_on_surface_fname);
		exit(0);
	}
	else{
		count = 0;
		while((fscanf(fr,"%lf%lf%lf",&tmp, &tmp, &T_on_surface[count]))!=EOF){
			count++;
		}
	}
	fclose(fr);

	if(strcmp(species_name,"xnu")!=0){
		printf("mu-surface name = %s\n", mu_on_surface_fname);
		fr = fopen(mu_on_surface_fname,"r");
		count = 0;
		while(fscanf(fr,"%lf %lf %lf", &tmp,&tmp,&mu_on_surface[count])!=EOF){
			count++;
		}
		fclose(fr);
	}

	/*
	 * Estimating the total xnu luminosity from the scattering surface.
	 */
	for(i=0; i<nrows_surf-1;i++){
		x0 = x[i]*1e5;
		z0 = z[i]*1e5;
		x1 = x[i+1]*1e5;
		z1 = z[i+1]*1e5;
		r0 = sqrt(x0*x0+z0*z0);// (1e5) is the conversion factor for km--->cm
		r1 = sqrt(x1*x1+z1*z1);// cm
		av_r = (r0+r1)/2.0;
		sintheta = (x0+x1)/(2.0*av_r);
		dl = sqrt((x1-x0)*(x1-x0)+(z1-z0)*(z1-z0));// cm
		dA = 2*M_PI*av_r*sintheta*dl; // Differential area (cm^2).
		av_dens = (density[i]+density[i+1])/2.0; //per cm^3.
		flux_xz[i] = (3.0/4.0)*c*av_dens*dA; // 1/s.
	}

	Nlumin_estimated_time = 0.0;
	for(i=0; i<nrows_surf-1;i++){
		Nlumin_estimated_time += flux_xz[i];
	}
	Nlumin_estimated_time = 2.0*Nlumin_estimated_time;
	
	/*
	 * Calculating the normalization factor for the x-neutrino number flux.
	 * normalization factor = NLuminosity provided in 'dd2_lumins_abs_vol_fname'/NLuminosity calculated from density on the surface.
	 */

	abs_dat abs_vol_dat  = get_L_abs_vol(time, species_name);
	Nlumin_provided_time = abs_vol_dat.NL_abs_vol;
	avE_provided_time    = abs_vol_dat.meanE_abs_vol;

	flux_renormalization_factor = Nlumin_provided_time/Nlumin_estimated_time;
	printf("Renormalization factor: %le\n",flux_renormalization_factor);

	/*
	Renormalizing the Density data on the surface and writting on to the file.
	*/
	strcpy(outfname,"../OutFiles/");
	strcat(outfname, species_name);
	strcat(outfname,"Density_percm3_on_surface_renormalized_1q_time.dat");

	fw = fopen(outfname,"w");


	for(i=0; i<nrows_surf; i++){
		if(strcmp(species_name, "xnu")==0){ //ormalizes the density on the surface inly for xnu.
			//printf("inside norm flux-sp_name: %s\n",species_name);
			density[i] = density[i]*flux_renormalization_factor;
		}
		fprintf(fw, "%lf\t%lf\t%le\n",x[i],z[i],density[i]);

	}
	fclose(fw);

	/*
	 * Recalculating the <Energy> of the x-neutrino on the scattering sutface and
	 * writting on to the file.
	 */

	outfname_avE_surf = (char *)calloc(150,sizeof(char));
	strcpy(outfname_avE_surf, "../OutFiles/");
	strcat(outfname_avE_surf, species_name);
	strcat(outfname_avE_surf, "avE_MeV_on_surface_renormalized_1q_time.dat");

	fw = fopen(outfname_avE_surf,"w");
	for(i=0; i<nrows_surf; i++){
		if(strcmp(species_name, "xnu")==0){
			avE_surface[i] = (1.0/2.0)*(3.15*T_on_surface[i]+avE_provided_time);
			//printf("%lf\t%lf\t%le\n",x[i],z[i],avE_surface[i]);
		}
		else{
			mubar = mu_on_surface[i]/T_on_surface[i];
			avE_surface[i] = T_on_surface[i]*(3.0*gsl_sf_fermi_dirac_int(3,mubar)/gsl_sf_fermi_dirac_2(mubar));
			//printf("%lf\t%lf\t%le\n",x[i],z[i],avE_surface[i]);
		}
		fprintf(fw, "%lf\t%lf\t%le\n",x[i],z[i],density[i]*avE_surface[i]);
	}
	fclose(fw);
	return(outfname_avE_surf);
}
