values_at_loc Manager(char *species_name, double time, double x_traj, double z_traj){

	char datfname[250];
	char *tmpfname;


	// ------------------ Intrpolating opticat depth data of the x-neutrino for the given time --------------
  //checked ---> NO SEG-FAULT !!!!

  strcpy(datfname, "tau");
  strcat(datfname, species_name);
  tmpfname = InterpolateData(datfname, time);
  char Taufname_time[200];
  strcpy(Taufname_time, tmpfname);

  // ------- Calculating the opticaldepth surface determined by tau = 2.0/3.0 for the x-neutrino corres ----
  //Checked ---> NO SEG-FAULT!!!

  tmpfname = TauSurface(Taufname_time, species_name);
  char Surfacefname_time[200];
  strcpy(Surfacefname_time, tmpfname);

  // ---------------------------------- Ordering the surface wrt theta ----------------------------------
  //Ckecked --> Removed a bug !!!!!!!

  tmpfname = OrderWrtTherta(Surfacefname_time,species_name,time);
  char SurfOrdered_time_1q[150];
  strcpy(SurfOrdered_time_1q, tmpfname);


  // ----------------------- Interpolating Temperature to the given time ---------------------------------
  //Checked ---> NO SEG-FAULT!!!

  strcpy(datfname, "T");
  tmpfname = InterpolateData(datfname, time);
  char Tfname_time[150];
  strcpy(Tfname_time, tmpfname);


  // -------------------- Interpolating chem.potential to the given time -------------------------------
  // xnu chem.pottential is assumed to be zero everywhere. So a null chem potential file is passed to
  //the subroutine.
  //Checked ---> NO SEG-FAULT!!!

  char mufname_time[150];
  if(strcmp(species_name,"xnu") != 0){
    strcpy(datfname, species_name);
    strcat(datfname, "mu");
    tmpfname = InterpolateData(datfname, time);
    strcpy(mufname_time, tmpfname);
  }
  else{
    strcpy(mufname_time, "\0");
  }

  // ----------------------------- Calculating the density of xneutrino assuming 0 chemical potential ---------
  //Checked ---> NO SEG-FAULT!!!

  tmpfname = DensityFD(Tfname_time, mufname_time);
  char Densityfname_time[150];
  strcpy(Densityfname_time, tmpfname);

  // ------------------------------ Interpolating nu species density on to the surface ---------------
  // Checked ---> No seg-fault..!!!!

  strcpy(datfname, species_name);
  strcat(datfname, "Density");
  tmpfname = ValueOnTheSurface(SurfOrdered_time_1q,Densityfname_time,datfname, time);
  char DensityOnSurffname_time[150];
  strcpy(DensityOnSurffname_time, tmpfname);

  // ------------------------------- Interpolating T density on to the tau surface -------------------
  // Checked ---> No seg-fault..!!!!

  strcpy(datfname, species_name);
  strcat(datfname, "T");
  tmpfname = ValueOnTheSurface(SurfOrdered_time_1q,Tfname_time,datfname, time);
  char TOnSurffname_time[150];
  strcpy(TOnSurffname_time, tmpfname);

  //----------------------------------- mu on surface --------------------------------------------------
  strcpy(datfname,species_name);
  strcat(datfname, "mu");
  tmpfname = ValueOnTheSurface(SurfOrdered_time_1q,mufname_time,datfname,time);
  char muOnSurffname_time[200];
  strcpy(muOnSurffname_time,tmpfname);

  // ------------------------------ Renormalizing the neutrino density -------------------------------- //
  // Checked ---> No seg-fault..!!!!

  tmpfname = NormalizeNFlux(species_name, DensityOnSurffname_time, TOnSurffname_time,muOnSurffname_time,time);
  char avE_surf_recal_fname_time[150];
  strcpy(avE_surf_recal_fname_time, tmpfname);
  free(tmpfname);	//created using malloc in NormalizeNFlux()
  char dens_on_surf_renorm_fname_time[150];
  strcpy(dens_on_surf_renorm_fname_time,outfname);

  // ------------------------------ Tracing the  number density at a point --------------------------
  //Checked ---> No seg-faults.
  tmpfname = profile_generator(species_name, dens_on_surf_renorm_fname_time, avE_surf_recal_fname_time,x_traj,z_traj);
  char traced_prof_fname_time[150];
  strcpy(traced_prof_fname_time,outfname);

  // ------------------------------------ ***************** ------------------------------------------

  values_at_loc values = Number_density_and_avEnergy(traced_prof_fname_time);
  printf("#%s @ given loc: %le, <E_%s> @ given loc: %le\n", species_name, values.n_nu, species_name,\
    values.av_e_nu);
  // ------------------------------------ ***************** ------------------------------------------
  
  return(values);

}
