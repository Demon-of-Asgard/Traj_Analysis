// Estimating values of given data on given surface (eg: Temparature values on neutrino surface)
char *ValueOnTheSurface(char *surffilename, char *datafilename, char* dataname, double time){

    strcpy(outfname,"../OutFiles/");
    strcat(outfname, dataname);
    strcat(outfname,"OnSurface_1q_time.dat");

    if(strcmp(datafilename, "\0")==0){
      printf("Datafilename is \\0. Returning!!\n");
      return(outfname);
    }
    
    else{
        FILE *fr,*fw;
        int i,j,k,count,surf_rows,surf_cols,N;
        dimf dim;
        double x1,x2,z1,z2, xcurrent,xnext,zcurrent,znext;
        double fxz,fx1z1,fx1z2,fx2z1,fx2z2,tmp,tmp2,x_loc,z_loc;
        double tmpdata[LNX*LNZ],data[LNX][LNZ];
        
        //printf("in_surf: %s\n",surffilename);
        //printf("in_data: %s\n",datafilename);
        dim  = shape_of_file(surffilename);
        surf_rows = dim.rows;
        surf_cols = dim.cols;
        //printf("#rows: %d,  #cols: %d in %s\n", surf_rows,surf_cols,surffilename);
        double x[surf_rows], z[surf_rows];
        double surfadata[surf_rows][3];
        count = 0;
        fr = fopen(surffilename,"r");//Reading the surface coordinates.
        while((fscanf(fr,"%lf%lf%lf",&x[count], &z[count], &tmp))!=EOF){
        count++;
        }
        fclose(fr);

        fr = fopen(datafilename,"r");//Reading the datafile to interpolate on to the surface coordinates.
        count = 0;
        while((fscanf(fr, "%lf%lf%lf", &tmp, &tmp,&tmpdata[count]))!=EOF){
        count++;
        }
        fclose(fr);

        //Converting above read 1D array of dim = LNX*LNZ to 2D array of dim = (LNX, LNZ).
        for(i=0; i<LNX; i++){
            for(j=0; j<LNZ; j++){
                N = i*LNZ+j;
                data[i][j] = tmpdata[N];
            }
        }

        fw = fopen(outfname,"w");

        //Interpolating data on to the surface.
        for(k=0;k<surf_rows;k++){
            x_loc=x[k];
            z_loc=z[k];

            for(i=0;i<LNX-1;i++){  
                xcurrent = (i-offset)*scale;
                xnext = (i+1-offset)*scale;

                if((xcurrent<=x_loc)&&(xnext>x_loc)){
                    for(j=0;j<LNZ-1;j++){
                        zcurrent = (j-offset)*scale;
                        znext = (j+1-offset)*scale;

                        if((zcurrent<=z_loc)&&(znext>z_loc)){ 
                            x1=xcurrent; x2=xnext; 
                            z1=zcurrent; z2=znext;
                            fx1z1=data[i][j];
                            fx1z2=data[i][j+1];
                            fx2z1=data[i+1][j]; 
                            fx2z2=data[i+1][j+1]; 
                            fxz = BilineatInterpolation(x_loc,z_loc,x1,z1,x2,z2,fx1z1,fx2z1,fx1z2,fx2z2);
                            fprintf(fw,"%lf\t%lf\t%lE\n",x_loc,z_loc,fxz);
                        }
                    }
                }
            }
        }
        fclose(fw);
        return(outfname);
    }
}