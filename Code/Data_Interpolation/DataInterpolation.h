/*
  Interpolate the given data, passed as an argument to 'dataname'
  to the time at 'time'.
  See 'Filenames.h' to see available datas.
  For e.g.,
  Consider the file name '../dd21352p5/dd21352p5tauxnu.dat',
  then the 'dataname' that should be passed on to the subroutine
  is tauxnu
  'Filenames.h' can be edited, if need be.
*/
char *InterpolateData(char *dataname, double time){
FILE *fr1, *fr2, *fw1;
char datafname1[150], datafname2[150];
char *fname1, *fname2;
int ix, iz, i,j,N;
double t1 = 0.0, t2 = 0.0;
double data1[LNX*LNZ];
double data2[LNX*LNZ];
double data_time[LNX*LNZ];
bool FLAG = true;
  if((time >=2.5)&&(time<5.0)){

    t1 = 2.5;
    t2 = 5.0;

    strcpy(datafname1,"../dd21352p5/dd21352p5");
    strcat(datafname1,dataname);
    strcat(datafname1,".dat");

    strcpy(datafname2, "../dd21355/dd21355");
    strcat(datafname2,dataname);
    strcat(datafname2, ".dat");

    fr1 = fopen(datafname1,"r");
    fr2 = fopen(datafname2,"r");

    if(fr1==NULL){
      FLAG = false;
      printf("%s does not exist.\n",datafname1);
    }
    if(fr2==NULL){
      FLAG = false;
      printf("%s does not exist.\n",datafname2);
    }
  }
  if((time >=5.0)&&(time<7.5)){

    t1 = 5.0;
    t2 = 7.5;

    strcpy(datafname1, "../dd21355/dd21355");
    strcat(datafname1,dataname);
    strcat(datafname1,".dat");

    strcpy(datafname2, "../dd21357p5/dd21357p5");
    strcat(datafname2,dataname);
    strcat(datafname2, ".dat");

    fr1 = fopen(datafname1,"r");
    fr2 = fopen(datafname2,"r");

    if(fr1==NULL){
      FLAG = false;
      printf("%s does not exist.\n",datafname1);
    }
    if(fr2==NULL){
      FLAG = false;
      printf("%s does not exist.\n",datafname2);
    }
  }
  if((time >=7.5)&&(time<=10.0)){

    t1 = 7.5;
    t2 = 10.0;

    strcpy(datafname1, "../dd21357p5/dd21357p5");
    strcat(datafname1,dataname);
    strcat(datafname1,".dat");

    strcpy(datafname2, "../dd213510/dd213510");
    strcat(datafname2,dataname);
    strcat(datafname2, ".dat");

    fr1 = fopen(datafname1,"r");
    fr2 = fopen(datafname2,"r");

    if(fr1==NULL){
      FLAG = false;
      printf("%s does not exist.\n",datafname1);
    }
    if(fr2==NULL){
      FLAG = false;
      printf("%s does not exist.\n",datafname2);
    }
  }
  if(FLAG == false){
    printf("File name(s) entered are incorrect. Check the subroutine 'double Interapolate_taux(double time)'\n");
    exit(0);
  }
  else{
    //printf("Interpolating %s data for the time:  %lf (ms)\n",dataname, time);
  }

  i = 0;
  while((fscanf(fr1,"%d%d%lf", &ix, &iz, &data1[i]))!=EOF){
    i++;
  }
  i = 0;
  while((fscanf(fr2,"%d%d%lf", &ix, &iz, &data2[i]))!=EOF){
    i++;
  }

  fclose(fr1);
  fclose(fr2);

  strcpy(outfname,"../OutFiles/");
  strcat(outfname, dataname);
  strcat(outfname,"_time.dat");

  fw1 = fopen(outfname, "w");
  for(i=0; i<LNX; i++){
    for(j=0; j<LNZ; j++){
      N = i*LNX+j;
      data_time[N] = NewtonsInterpolation(time, t1, t2, data1[N], data2[N]);
      fprintf(fw1,"%d\t%d\t%le\n", i+1, j+1,data_time[N]);
    }
  }
  fclose(fw1);
  return(outfname);
}
