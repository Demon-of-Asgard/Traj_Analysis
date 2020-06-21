/*
 * Builtin header files.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_fermi_dirac.h>

/*
 *Order of inclusion of the user defined header files is very important.
 */
#include "../Global/Global.h"
#include "../Filenames/Filenames.h"
#include "../Subsidiary_Subroutines/shape_of_file.h"
#include "../Subsidiary_Subroutines/Interpolation_subroutines.h"
#include "../Data_Interpolation/DataInterpolation.h"
#include "../Optical_Depth_Surface/OpticaldepthSurface.h"
#include "../Subsidiary_Subroutines/Order_wrtTheta.h"
#include "../Density/Density.h"
#include "../Value_On_the_Surface/ValueOnTheSurface.h"
#include "../Get_L_abs_vol/Get_L_abs_vol.h"
#include "../Norm_Flux/NormNFlux.h"
#include "../Tracer/Tracer.h"
#include "../Integrate_Spectrum/Integrate_spectrum.h"
#include "../Binterpole/binterpole.h"
#include "../Manager/Manager.h"
#include "../Master/Master.h"
//#include "rmfile.h"