//==========================================================================================================================================

double NewtonsInterpolation(double x, double x1, double x2, double f1, double f2){

  double slope;
  double f;
  slope = (f2-f1)/(x2-x1);
  f = f1 + slope*(x-x1);
  return(f);
}

//==========================================================================================================================================

double InverseNewtonsInterpolation(double x1, double x2, double f1, double f2, double f){

  double slope = 0.;
  double x = 0.;
  slope = (f2-f1)/(x2-x1);
  x = x1+((f-f1)/slope);
  //x = (f - (f1 * x2 - f2 * x1) / (x2 - x1)) / ((f2 - f1) / (x2 - x1));
  return(x);
}

//==========================================================================================================================================


double BilineatInterpolation(double x,double z,double x1,double z1,double x2,double z2,double fx1z1,double fx2z1,double fx1z2,double fx2z2){

  double fxz=(1/((x2-x1)*(z2-z1)))*((fx1z1*(x2-x)*(z2-z))-(fx2z1*(x1-x)*(z2-z))-(fx1z2*(x2-x)*(z1-z))+(fx2z2*(x1-x)*(z1-z)));
  return(fxz);
}

//==========================================================================================================================================


double LogLinInterpolation(double x, double x1, double x2, double f1, double f2){
    
    double logf = 0.0;
    double slope = (log(f2)-log(f1))/(x2-x1);

    if(x>=50.0 && x<=100.0){ 
        logf = log(f1) + slope*(x-x1);
    }
    
    else if(x>100){
        logf = log(f2) + slope*(x-x2);
    }

    else if(x<50.0){
        logf = log(f1) - slope*(x-x1);
    }
    double f = exp(logf);
    return(f);
}

//==========================================================================================================================================

double LinInverseR2Interpolation(double x, double x1, double x2, double f1, double f2){

  double a = ((f2*x2*x2)-(f1*x1*x1))/((x2*x2)-(x1*x1));
  double b = (f2-f1)*(x1*x1*x2*x2)/((x1*x1)-(x2*x2));
  double f = a + b/(x*x);
  return(f);
}

//==========================================================================================================================================
