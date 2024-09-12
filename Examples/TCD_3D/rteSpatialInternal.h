
void ExampleFile()
{
  
  OutPut("Example: rteSpatial.h" << endl);  
}

// ========================================================================
// definitions for the temperature
// ========================================================================

void Exact_RTE( double *X, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double k = TDatabase::ParamDB->P0, x, y, z, s1, s2, s3;

  x = X[0];
  y = X[1]; 
  z = X[2];
  s1 = X[3];
  s2 = X[4];
  s3 = X[5];

  values[0] = (exp(-k*t))*s1*s2*s3*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
}

void Exact( double x, double y, double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double k = TDatabase::ParamDB->P0;
  
  double s1 = TDatabase::ParamDB->P11;
  double s2 = TDatabase::ParamDB->P12;
  double s3 = TDatabase::ParamDB->P13;

  values[0] = (exp(-k*t))*s1*s2*s3*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[1] = Pi*(exp(-k*t))*s1*s2*s3*cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[2] = Pi*(exp(-k*t))*s1*s2*s3*sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[3] = Pi*(exp(-k*t))*s1*s2*s3*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*s1*s2*s3*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
}
 

void BoundCondition(int dummy,double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int dummy,double x, double y, double z, double &value)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double k = 0.1;
  
  value = (exp(-k*t))*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);  
}

void GetKernel(int N_Inputs, double *Inn, double *Out) 
{

 double c1 = 1./(4.*Pi);
 Out[0] =  c1;
//  Out[0] = c1*(1.+ (Inn[0]*Inn[3] + Inn[1]*Inn[4] + Inn[2]*Inn[5]));

}

void GetRhs(int N_Inputs, double *Inn, double *Out) 
{
  double s1 = Inn[0];
  double s2 = Inn[1];
  double s3 = Inn[2];
  double x= Inn[6];
  double y= Inn[7];
  double z= Inn[8];
  double alpha = TDatabase::ParamDB->P0;  
  double sigma_a = TDatabase::ParamDB->P1;
  double sigma_s = TDatabase::ParamDB->P2;
  double spx = sin(Pi*x);
  double spy = sin(Pi*y);
  double spz = sin(Pi*z);
  double expt = exp(-alpha*TDatabase::TimeDB->CURRENTTIME);

 Out[0] = (sigma_a - sigma_s)*s1*s2*s3*expt*spx*spy*spz; //  

}
// ========================================================================
// BilinearCoeffs for RTE in x-direction
// ========================================================================
void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, alpha, sigma_a, sigma_s, spx, spy, spz;                                  // *param;
  double x, y, z, c, a[3],s[3], h;
  double expt;
  
  alpha = TDatabase::ParamDB->P0;
  sigma_a = TDatabase::ParamDB->P1;
  sigma_s = TDatabase::ParamDB->P2;
  expt = exp(-alpha*TDatabase::TimeDB->CURRENTTIME);
  double diffusion, PE_NR = TDatabase::ParamDB->PE_NR;

  if(PE_NR==0.)
   { diffusion = 0;}
  else
   {diffusion = 1./PE_NR;}

  alpha = -alpha + 3.*Pi*Pi*diffusion;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];
    z = Z[i];

    spx = sin(Pi*x);
    spy = sin(Pi*y);
    spz = sin(Pi*z);
   
    // diffusion
    coeff[0] = diffusion;
    // reaction
    coeff[4] = 0.;

     // rhs
    coeff[1] = Pi*cos(Pi*x)*spy*spz*expt;  // s1 part
    coeff[2] = Pi*spx*cos(Pi*y)*spz*expt;  // s2 part
    coeff[3] = Pi*spx*spy*cos(Pi*z)*expt;  // s3 part
    coeff[5] = alpha*expt*spx*spy*spz;  // f
  }
}

void SystemRhs(int N_Active, double *param, double *B, double *p1_null, double **RHSs, double **p2_null)
{
 double s1, s2, s3, smult, tau;
  
  s1 = param[0];
  s2 = param[1];  
  s3 = param[2];
  smult = s1*s2*s3;
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 

  Daxpy(N_Active, smult*tau*TDatabase::TimeDB->THETA4,  RHSs[0], B);    
  Daxpy(N_Active, s1*smult*tau*TDatabase::TimeDB->THETA4, RHSs[1], B);   
  Daxpy(N_Active, s2*smult*tau*TDatabase::TimeDB->THETA4,  RHSs[2], B);   
  Daxpy(N_Active, s3*smult*tau*TDatabase::TimeDB->THETA4,  RHSs[3], B);  

  //supg 
   if (TDatabase::ParamDB->DISCTYPE==SUPG)
    {
     Daxpy(N_Active, s1*smult*tau,  RHSs[4], B);   
     Daxpy(N_Active, s2*smult*tau,  RHSs[5], B);   
     Daxpy(N_Active, s3*smult*tau,  RHSs[6], B);  

     Daxpy(N_Active, s1*s1*smult*tau,  RHSs[7], B);   
     Daxpy(N_Active, s2*s2*smult*tau,  RHSs[8], B);   
     Daxpy(N_Active, s3*s3*smult*tau,  RHSs[9], B);  

     Daxpy(N_Active, s1*s2*smult*tau,  RHSs[10], B);   
     Daxpy(N_Active, s1*s3*smult*tau,  RHSs[11], B);   
     Daxpy(N_Active, s2*s3*smult*tau,  RHSs[12], B);  
    }

}


// kind of boundary condition (for FE space needed)
void SurfBoundCondition(int BdComp, double t, BoundCond &cond)
{
//   cond = DIRICHLET;
//  each edge of all triangles on a 3D surface is an interior edge
  cout << "Error!! each edge of all triangles on a 3D surface is an interior edge"<<endl;
  cout << "Error!! Check the Example file "<<endl;
  exit(0);

}

// value of boundary condition
void SurfBoundValue(int BdComp, double Param, double &value)
{
//   value = 0;
//  each edge of all triangles on a 3D surface is an interior edge
  cout << "Error!! each edge of all surface triangles on a 3D surface is an interior edge"<<endl;
  cout << "Error!! Check the Example file "<<endl;
  exit(0);
}

void Initial_All(int N_Coord, double* in, double* out)
{
  
}


void InitialS(double x, double y, double *values)
{
   double h, r;

 double t = TDatabase::TimeDB->CURRENTTIME;
//  if(fabs(x)<1e-5  && fabs(y)<1e-5 )
//  cout << "InitialS : x : " << x << " y : " << y <<endl;

  values[0] = exp(-6.*t)*x*y;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  cout << "Error!! Check the Example file "<<endl;
  exit(0);

}
void ExactS(double x, double y, double *values)
{


  cout << "Error!! Check the Example file "<<endl;
  exit(0);

 double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = exp(-6.*t)*x*y;
  values[1] = exp(-6.*t)*y;
  values[2] = exp(-6.*t)*x;
  values[3] = 0;
  values[4] = 0;

}

void InitialSall(double x, double y,  double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = x*y*exp(-6.*t);
}

void ExactSall_oldStep(double x, double y,  double z, double *values)
{
 double t = TDatabase::ParamDB->P14; // old time step should be storen in main prgram

  values[0] = x*y*exp(-6.*t);
  values[1] = exp(-6.*t)*y;
  values[2] = exp(-6.*t)*x;
  values[3] = 0;
  values[4] = 0;
}


void ExactSall(double x, double y,  double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = x*y*exp(-6.*t);
  values[1] = exp(-6.*t)*x;
  values[2] = exp(-6.*t)*y;
  values[3] = 0;
  values[4] = 0;
}

void SurfAllCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1.;
  int i;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

  }

  cout << "Error!! Check the Example file "<<endl;
  exit(0);

}


void SurfCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double r2;
  static double eps = 1./TDatabase::ParamDB->RE_NR;

  for(i=0;i<n_points;i++)
   {
    coeff = coeffs[i];
    coeff[0] = eps;   // eps

    if(TDatabase::ParamDB->FR_NR == 0)
       coeff[1] = 0;
    else
       coeff[1] = TDatabase::ParamDB->FR_NR;

   }
  cout << "Error!! Check the Example file "<<endl;
  exit(0);
}


// For Siminhale 
void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
//   cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
}

void BoundValue_LMin(double x, double y, double z, double *values)
 {
   
   cout<< "Nucleation implemented separately " << endl;
   exit(0);

   //nucleation
   // no nucleation, so f = 0 at L=0
   
//    double t;   
//  if((x==0) && y>(1./3) &&  y<(2./3))
//   {
//   if((x==0) && y>(1./3) &&  y<(1./3. + 1./24.))
//    {
//     t = -6. + 12.*(y - 1./3.)*24.*(1. + 1./3. + 1./24. - y);
//     t = 1./(1. + exp(-t));
// //              cout << y << " t " << t << endl;
//     values[0] = t;
//     values[1] = 0.;
//     values[2] = 0.;    
//    }
//   else if((x==0) && y>(2./3 - 1./24.)  &&  y<(2./3))
//    {
//      t = 6. + 12.*(2./3. - 1./24. - y)*24.*(1. + 2./3. -y);
//      t = 1./(1. + exp(-t));
// //          cout << y << " t " << t << endl;
//     values[0] = t;
//     values[1] = 0.;
//     values[2] = 0.;    
//    }  
//   else
//    {
//     values[0] = 1.;
//     values[1] = 0.;
//     values[2] = 0.;
//    }
//   }
//  else
//   {
//     values[0] = 0.;
//     values[1] = 0.;
//     values[2] = 0.;
//   } 

//   double L_min = TDatabase::ParamDB->UREA_D_P_0;
//   double L_max = TDatabase::ParamDB->UREA_D_P_MAX;
//   double f_infty = TDatabase::ParamDB->UREA_f_infty;
//   double eps = 1e-8, val, a_min, g;
//   double f_in = TDatabase::ParamDB->UREA_f_infty;
// 

 }


void BoundValue_LMax(double x, double y,  double z,  double *values)
 {
    values[0] = 0.;
    values[1] = 0.;
    values[2] = 0.; 
 }
