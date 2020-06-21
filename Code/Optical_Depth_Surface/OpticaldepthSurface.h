char *TauSurface(char *infname,char *type){
  FILE *fr, *fw;
  strcpy(outfname,"../OutFiles/");
  strcat(outfname,type);
  strcat(outfname,"tau_surface_time.dat");
  int i,j,N, tmpi;
  double SurfValue = 2.0/3.0;
  double tau_tmp[LNX*LNZ];
  double tau[LNX][LNZ];
  double x1, x2,x;
  double data1, data2;

  //Reading the tau file.
  fr = fopen(infname, "r");
  if(fr==NULL){
    printf("Unable to open %s, breaking :(\n",infname);
    exit(0);
  }
  i = 0;
  while((fscanf(fr,"%d%d%lf",&tmpi, &tmpi, &tau_tmp[i]))!=EOF){
    i++;
  }
  fclose(fr);
  //Changing 1D array to 2D array.
  for(i=0;i<LNX;i++){
    for(j=0;j<LNZ;j++){
      N = i*LNZ+j;
      tau[i][j] = tau_tmp[N];
    }
  }

  // Calculating optical depth surface.
  int count = 0;
  fw = fopen(outfname,"w");
  if(fr==NULL){
    printf("Unable to open %s, breaking :(\n",outfname);
    exit(0);
  }
  for(i=1;i<LNX-1;i++){
    for(j=1;j<LNZ-1;j++){

      if(((tau[i][j]>=SurfValue)&&(tau[i][j+1]<=SurfValue))||((tau[i][j]<=SurfValue)&&(tau[i][j+1]>=SurfValue))){
        x1 = (double)(j+1);
        x2 = (double)(j+2);
        data1 = tau[i][j];
        data2 = tau[i][j+1];
        x = InverseNewtonsInterpolation(x1,x2,data1,data2,SurfValue);
        count ++;
        fprintf(fw,"%2.5lf\t %2.5lf\t %le\n",(i+1-145.0)*scale,(x-145.0)*scale,SurfValue);
      }

      else if(((tau[i][j]>=SurfValue)&&(tau[i+1][j]<=SurfValue))||((tau[i][j]<=SurfValue)&&(tau[i+1][j]>=SurfValue))){
        x1 = (double)(i+1);
        x2 = (double)(i+2);
        data1 = tau[i][j];
        data2 = tau[i+1][j];
        x = InverseNewtonsInterpolation(x1,x2,data1,data2,SurfValue);
        count ++;
        fprintf(fw,"%2.5lf\t %2.4lf\t %le\n",(x-145.0)*scale,(j+1-145.0)*scale,SurfValue);
      }

      else{
        continue;
      }
    }
  }
  fclose(fw);
  return(outfname);
}
