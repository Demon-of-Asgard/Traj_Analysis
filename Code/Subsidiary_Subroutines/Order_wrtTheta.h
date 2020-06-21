char *OrderWrtTherta(char *infilename, char* species_name, double time){
  //printf("in file name : %s \n", infilename);

  strcpy(outfname,"../OutFiles/");
  strcat(outfname,species_name);
  strcat(outfname,"SurfOrdered_time_1q.dat");
  FILE *fr, *fw;
  int i,j, nrows = 0, ncols = 0;
  double x_temp,z_temp,data_temp,theta_temp;

  // counting #columns and rows.
  dimf dim = shape_of_file(infilename);
  nrows = dim.rows;
  ncols = dim.cols;
  //printf("#rows: %d, #columns: %d in %s.\n", nrows, ncols, infilename);
  double xx,zz, rr, angle, x[nrows], z[nrows], r[nrows], theta[nrows], data[nrows];

  /*----------------- Calculating the angle -----------------*/
  fr = fopen(infilename, "r");
  i = 0;
  while((fscanf(fr, "%lf%lf%lf",&x[i], &z[i], &data[i]))!=EOF){
    xx = fabs(x[i]);
    zz = fabs(z[i]);
    rr = sqrt(xx*xx+zz*zz);
    angle = acos(zz/rr);
    r[i] = rr;

    if(x[i]==0. && z[i]>0.)
      {theta[i] = 0.0;}
    if(x[i]>0. && z[i]>0)
      {theta[i] = angle;}
    if(x[i]>0  && z[i] == 0)
      {theta[i] = M_PI/2;}
    if(x[i]>0. && z[i] <0)
      {theta[i] = (M_PI)-angle;}
    if(x[i]==0. && z[i]<0)
      {theta[i] = M_PI;}
    if(x[i]<0. && z[i]<0.)
      {theta[i] = M_PI+angle;}
    if(x[i]<0 && z[i] == 0)
      {theta[i] = 3.*M_PI/2.;}
    if(x[i]<0. && z[i]>0.)
      {theta[i] = (2.*M_PI)-angle;}
    i++;
  }
  fclose(fr);

  /*----------------- Ordering w.r.t angle -----------------*/
  for(i=0; i<nrows; i++){
    for(j=0; j<nrows; j++){
      if(theta[j]>theta[i]){
        theta_temp = theta[i];
        theta[i] = theta[j];
        theta[j] = theta_temp;

        x_temp = x[i];
        x[i] = x[j];
        x[j] = x_temp;

        z_temp = z[i];
        z[i] = z[j];
        z[j] = z_temp;

        data_temp = data[i];
        data[i] = data[j];
        data[j] = data_temp;
      }
    }
  }

  /*----------------- Writting to the file. -----------------*/
  fw = fopen(outfname, "w");
  for(i = 0; i<nrows; i++){
    if((theta[i]>=0.)&&(theta[i]<=M_PI/2.0)){
      if(sqrt(x[i]*x[i]+z[i]*z[i])<75.0){
        fprintf(fw,"%lf\t%lf\t%lE\n", x[i], z[i],data[i]);
      }
    }
  }
  fclose(fw);
  return(outfname);
}
