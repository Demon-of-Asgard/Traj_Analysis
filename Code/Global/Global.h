//Reaytrace properties
#define NW 1   // number of omega[:=1/(neutrino energy)] grids
#define NU 400 // number of costheta grids u=cos\theta
#define NP 400 // # of phi grids
#define NX 1
#define NZ 1

//bin proparties
#define nbins 18

//Grid Properties
int LNX = 289;
int LNZ = 289;
double offset = 145.0;
double scale = 0.739;
double trajtime_offset = 0.1700E+02; //trajtime_offset; //(3441*0.004926);


double c = 3.0e10; // (cm/s)

char outfname[150];

// dimension of files ---> dimf.
typedef struct{
	unsigned int rows;
	unsigned int cols;
}dimf;

//result of xnu ray tracing
typedef struct{
	double n_nu;
	double av_e_nu;
}values_at_loc;

//Required to interpolate the bin data.
typedef struct{
	double Nnu_rt;
	double Nanu_rt;
	double avEnu_rt;
	double avEanu_rt;
}bin_polated; // Extrapolated using bin data.

//Required by Get_L_abs_vol.h
typedef struct{
	double NL_abs_vol;
	double meanE_abs_vol;
}abs_dat;
