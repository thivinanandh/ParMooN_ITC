// Brennstoffzelle

void ExampleFile()
{
  OutPut("Example: Brennstoffzelle.h" << endl) ;
}

// ========================================================================
// initial condition
// ========================================================================
void InitialU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
   //if ((x>=40) && (x<=50) && (y>=30) && (y<=40) && (fabs(z-5)<=1e-6))
   if ((((x-48)*(x-48)+(y-20)*(y-20)-16 ) <=1e-6) && (fabs(z-5)<=1e-6))
   {
      cond = NEUMANN;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
   }
   else
      cond  = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
   //if ((x>=10) && (x<=20) && (y>=100) && (y<=110) && (fabs(z-5)<=1e-6))
   if ((((x-17)*(x-17)+(y-120)*(y-120)-16 ) <=1e-6) && (fabs(z-5)<=1e-6))
   {
      value = -1;
   }
   else
   {
      value = 0;
   }
   // if ((((x-48)*(x-48)+(y-20)*(y-20)-16 ) <=1e-6) && (fabs(z-5)<=1e-6))
   //{
   //   value = -1;
   // }

}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
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
}
