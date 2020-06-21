char *DensityFD(char *Tfname, char *Mufname){

  bool flag1 = false, flag2 = false;
  int i, j,N,count;
  double T[LNX][LNZ], mu[LNX][LNZ],densityMeV3,densitypercm3, x[LNX*LNZ], z[LNX*LNZ], data[LNX*LNZ],mubar;

  FILE *fr, *fw;
  /*----------------- Reading T -----------------*/
  if (strcmp(Tfname, "\0")==0){
    printf("Null T file.\n");
  }
  else{
    flag1 = true;
    fr = fopen(Tfname, "r");
    count =0;
    while((fscanf(fr,"%lf%lf%lf",&x[count],&z[count],&data[count]))!=EOF){
      count++;
    }
    fclose(fr);
    for(i=0;i<LNX;i++){
      for(j=0;j<LNZ;j++){
        N = i*LNZ+j;
        T[i][j] = data[N];
      }
    }
  }

  /*----------------- Reading mu -----------------*/
  if(strcmp(Mufname, "\0")== 0){
    //printf("Null mu file.\n");
    for(i=0;i<LNX;i++){
      for(j=0;j<LNZ;j++){
        N = i*LNZ+j;
        mu[i][j] = 0.0;
      }
    }
  }
  else{
    flag2 = true;
    fr = fopen(Mufname,"r");
    count = 0;
    while((fscanf(fr,"%lf%lf%lf",&x[count],&z[count],&data[count]))!=EOF){
      count++;
    }
    fclose(fr);
    for(i=0;i<LNX;i++){
      for(j=0;j<LNZ;j++){
        N = i*LNZ+j;
        mu[i][j] = data[N];
      }
    }
  }

  /*-------- Calculating density (cm)^-3 and Writting to the density file --------*/
  if (flag1 || flag2){
    strcpy(outfname, "../OutFiles/Densitypercm3_time.dat");
    fw = fopen(outfname, "w");
    if(fw==NULL){
      printf("Unable to open %s. Exiting :(\n", outfname);
    }
    for(i=0; i<LNX; i++){
      for(j=0; j<LNZ; j++){
        mubar = mu[i][j]/T[i][j];
        densityMeV3 = (T[i][j] * T[i][j] * T[i][j])*(2.0/ (2.0*M_PI*M_PI))*gsl_sf_fermi_dirac_2(mubar);
        densitypercm3 = (1.3*1E32)*densityMeV3;
        fprintf(fw, "%d\t%d\t%le\n", i+1,j+1,densitypercm3);
      }
    }
    fclose(fw);
  }
  else{
    printf("Entered NULL files. Unable to calculate density.\n");
  }
  return(outfname);
}