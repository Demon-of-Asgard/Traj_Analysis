#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#define nbins 18
#define datperbin 4

int switchfile(int filer){
	printf("\nSwitching files\n");
	printf("filer = %d\n", filer);
	if (filer==50){
		return (100);
	}
	if (filer ==100){
		return(50);
	}
}

int switchspecies(int species){
	printf("\nSwitching species\n");
	printf("species = %d\n", species);
	if (species==1){
		return (2);
	}
	if (species ==2){
		return(1);
	}
}

int main(){
	FILE *fr,*fw, *fwL50, *fwL100, *fwE50, *fwE100;
	FILE *fwLnu50,*fwLanu50,*fwLnu100,*fwLanu100,*fwEnu50,*fwEanu50,*fwEnu100,*fwEanu100;
	double dataperrow[datperbin][nbins], tmpdat, time;
	int count = 0, datkind =1; //datakind ---> 1: \nu_e, 2: \bar\nu_e 
	int bin = 0;
	int Ndat=0;
	bool atEOF = false;
	int file = 50;
	int species = 1;
//Opening the o/p files
	fwL50 = fopen("lumin50.dat", "w");
	fwL100 = fopen("lumin100.dat","w");
	fwE50 = fopen("avE50.dat", "w");
	fwE100 = fopen("avE100.dat","w");

	fwLnu50 = fopen("luminnu50.dat", "w");
	fwLnu100 = fopen("luminnu100.dat","w");
	fwLanu50 = fopen("luminanu50.dat", "w");
	fwLanu100 = fopen("luminanu100.dat","w");

	fwEnu50 = fopen("avEnu50.dat", "w");
	fwEanu50 = fopen("avEanu50.dat", "w");
	fwEnu100 = fopen("avEnu100.dat","w");
	fwEanu100 = fopen("avEanu100.dat","w");
// opening the i/p file
	fr = fopen("../dd2_135_135_nu16420/neutrino_properties_bins.dat", "r");
	if(!fr){
		printf("Null pointer\n");
	}
//Reading the i/p file
	else{
		fscanf(fr, "%lE%lE", &tmpdat, &tmpdat);
		while((fscanf(fr, "%lE", &tmpdat))!=EOF)
		{
			count++;
		}
		printf("Length of the dat file: %d\n", count);
	}
	fclose(fr);

	count = 0;
	fr = fopen("../dd2_135_135_nu16420/neutrino_properties_bins.dat", "r");
	if(!fr){
		printf("Null pointer\n");
		exit(0);
	}
	else{
		fscanf(fr, "%lE%lE", &tmpdat, &tmpdat);// first two dat ------->#bins & angular width per bin
//Reading and entering the staring time of the i/p file  into the o/p files
		if(fscanf(fr, "%lE", &time) != EOF){
			fprintf(fwL50, "%lE\t", time);
			fprintf(fwL100, "%lE\t", time);
			fprintf(fwE50, "%lE\t", time);
			fprintf(fwE100, "%lE\t", time);

			fprintf(fwLnu50, "%lE\t", time);
			fprintf(fwLnu100, "%lE\t", time);
			fprintf(fwEnu50, "%lE\t", time);
			fprintf(fwEnu100, "%lE\t", time);

			fprintf(fwLanu50, "%lE\t", time);
			fprintf(fwLanu100, "%lE\t", time);
			fprintf(fwEanu50, "%lE\t", time);
			fprintf(fwEanu100, "%lE\t", time);
		}
//Reading and witting the rest of the file
		while(fscanf(fr, "%lE", &tmpdat) != EOF){
			count++;

			if (datkind == 1){ // datakind = 1 for luminosity and 2 for av. Energy.
				if (file == 50){ // file == 50 for bins at radius 50 and 100 for bins at radius 100.
					fprintf(fwL50, "%lE\t", tmpdat);
					if(species == 1){ // species =1 for nu and 2 for anu.
						fprintf(fwLnu50, "%lE\t", tmpdat);
					}
					if(species == 2){
						fprintf(fwLanu50, "%lE\t", tmpdat);
					}
				}
				if (file == 100){
					fprintf(fwL100, "%lE\t", tmpdat);
					if(species == 1){
						fprintf(fwLnu100, "%lE\t", tmpdat);
					}
					if(species == 2){
						fprintf(fwLanu100, "%lE\t", tmpdat);
					}
				}
			}

			if (datkind == 2){
				if (file == 50){
					fprintf(fwE50, "%lE\t", tmpdat);
					if(species == 1){
						fprintf(fwEnu50, "%lE\t", tmpdat);
					}
					if(species == 2){
						fprintf(fwEanu50, "%lE\t", tmpdat);
					}
				}
				if (file == 100){
					fprintf(fwE100, "%lE\t", tmpdat);
					if(species == 1){
						fprintf(fwEnu100, "%lE\t", tmpdat);
					}
					if(species == 2){
						fprintf(fwEanu100, "%lE\t", tmpdat);
					}
				}
			}

			if((count%36) == 0){ // switching to the radius 100 bins after reading the 36 datas (#bins*#species) and back after reading next 36.
				file =  switchfile(file);
				printf("file: %d\n", file);
			}

			if (count%18 == 0){ // Switchinfg the species after reading each 18 data.
				species = switchspecies(species);
			}
			if (count == 72){ // switching the data kind after reading every 72(#bins*#species*#positions (r=50 and r=100)).
				datkind = 2;
				printf("%d\t%d\n",count, datkind);
			}
			if(count == 144){// Resetting the datkind afterreading the first row
				datkind=1;
				printf("%d\t%d\n",count, datkind);
			}

			if(count == 144){//moving to the next row after the completion of each row.
				fprintf(fwL50, "\n");
				fprintf(fwL100, "\n");
				fprintf(fwE50, "\n");
				fprintf(fwE100, "\n");

				fprintf(fwLnu50, "\n");
				fprintf(fwLnu100, "\n");
				fprintf(fwEnu50, "\n");
				fprintf(fwEnu100, "\n");

				fprintf(fwLanu50, "\n");
				fprintf(fwLanu100, "\n");
				fprintf(fwEanu50, "\n");
				fprintf(fwEanu100, "\n");

				if(fscanf(fr, "%lE", &time) != EOF){ // reading the first entry of eact row from i/p file and writting them to the o/p files

					fprintf(fwL50, "%lE\t", time);
					fprintf(fwL100, "%lE\t", time);
					fprintf(fwE50, "%lE\t", time);
					fprintf(fwE100, "%lE\t", time);

					fprintf(fwLnu50, "%lE\t", time);
					fprintf(fwLnu100, "%lE\t", time);
					fprintf(fwEnu50, "%lE\t", time);
					fprintf(fwEnu100, "%lE\t", time);

					fprintf(fwLanu50, "%lE\t", time);
					fprintf(fwLanu100, "%lE\t", time);
					fprintf(fwEanu50, "%lE\t", time);
					fprintf(fwEanu100, "%lE\t", time);

					printf("%d\t%lE\n", count, time);
					count =0;
				}

			}

		}

	}
fclose(fwL50);
fclose(fwL100);
fclose(fwE50);
fclose(fwE100);

fclose(fwLnu50);
fclose(fwLnu100);
fclose(fwEnu50);
fclose(fwEnu100);
fclose(fwLanu50);
fclose(fwLanu100);
fclose(fwEanu50);
fclose(fwEanu100);
return (0);
}