// channel with circular cross section



#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>
#include <IsoBoundFace.h>
#include <MacroCell.h>
#include <BdSphere.h>
#include <tetgen.h>


void ExampleFile()
{
  OutPut("Example: CircularChannel.h" << endl) ;
  
//   #define __Cylinder__   
}


// For Bent Pipe cylinder

// 1000 - Inlet 
// 1001 - wall
// 1002 - Outlet

// Inlet facing in yZ plane and velocity inlet is in x direction
// The pipe starting at z = 0 and bending towards z = 10 ( positive z direction 
// So the "g" value of gravity is 9.81 m/s^2 (positive z direction)

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
 
  values[0] = 1-(x*x+y*y);
  values[1] = -2*x;
  values[2] = -2*y;
  values[3] = 0;
  values[4] = -4;
}

void ExactP(double x, double y,  double z, double *values)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;

  values[0] = 4*eps*(10-z);
  values[1] = 0;
  values[2] = 0;
  values[3] = -4*eps;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
 // (1000) Inlet - dirichlet, (1001) wall - dirichlet, (1002) Outlet - neumann
    if(CompID == 0)
    {
        cond = DIRICHLET;
    }
    else if(CompID == 1)
    {
        cond = DIRICHLET;
    }
    else if(CompID == 2)
    {
        cond = NEUMANN;
    }
    else
    {
        Error("Unknown Boundary component");
    }
}


// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
    // INlet : value is 1, everywhere else 0
    if(CompID == 0)
    {
        // if the point is less than 0.9 times of radius of inlet, then set the value to 1
        // the inlet is on yz plane, centered at 0,0 with radius 0.5
        double radius = 0.5;
        double dist = sqrt(y*y + z*z);
        if( dist < 0.9*radius)
            value = 1;
        else
            value = 0;
    }
    else
        value = 0;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
   value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
      
    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = -10.19367;
    coeff[3] = 0;  
    }
}

void ReadMeditMesh(char *SMESH, tetgenio &In)
{  
  
//    exit(0);    
} // ReadMeditMesh


void TetrameshGen(TDomain *&Domain)
{
 // Removed , as its not used

}



