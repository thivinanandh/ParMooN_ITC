#include<NodalFunctional1D.h>
#include<RefTrans1D.h>
#include<LineAffin.h>
#include<MacroCell.h>
#include<HexaAffin.h>
#include<QuadAffin.h>
#include<TetraAffin.h>


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


void InitialValue(int N_Inputs, double *Inn, double *Out) 
{
  // Inputs X, Y, Z, l1
  // outputs U
  double x = Inn[0];
  double y = Inn[1];
  double z = Inn[2];
  double l1 = Inn[3];
  Out[0] = l1;
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


//// ------------- SIMINHALE THIVIN -----------------////

void BoundCondition_Velocity(int BdComp,  double x, double y, double z,  BoundCond &cond)
{
  	if (BdComp == 0 || BdComp == 2)
		cond = DIRICHLET;

	else
	{
		cond = NEUMANN;
	}
}

// Boundary Condition for the internal System 
void BoundCondition_LminLMax_L0(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;

  // cond_Lmin = DIRICHLET;
  // cond_Lmax = DIRICHLET;
}


// For Siminhale 
 void BoundCondition_LminLMax_L1(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
  // cond_Lmin = DIRICHLET;
  // cond_Lmax = DIRICHLET;  
}


void BoundValue_LMin(double x, double y, double z, double *values)
 {
   
   cout<< "Nucleation implemented separately " << endl;
   exit(0);
 }


void BoundValue_LMax(double x, double y,  double z,  double *values)
 {
    values[0] = 0.;
    values[1] = 0.;
    values[2] = 0.; 
 }



// Nucleation Kernel
void GrowthAndB_Nuc_L0(int N_Inputs, double *Inn, double *Out) 
{
  double B_Nuc;
  double R_0  = TDatabase::ParamDB->P6; // Nucleation factor
  // double R_0  =3.35*1.01663385882;
  double f1  = TDatabase::ParamDB->P7; // interactive index
  // double bsigma  = TDatabase::ParamDB->P8;     
  // double sigma_c = TDatabase::ParamDB->P9;  
  double vsymp = TDatabase::ParamDB->P1;
  // double health  = TDatabase::ParamDB->P10;    
  // double bhealth = TDatabase::ParamDB->P11;  
  // double health_c = TDatabase::ParamDB->P12;  
  double sd  = TDatabase::ParamDB->P13; 
  double a4= 4.7, b4=0.28, c4_sqr=0.12*0.12;      
  // double t=TDatabase::TimeDB->CURRENTTIME ; 

  if(N_Inputs!=8)
  {
    printf("N_Inputs != 8,  Get_GrowthAndB_Nuc: %d \n",N_Inputs); 
#ifdef _MPI
     MPI_Finalize();
#endif       
    exit(0);        
  }

  // double x = Inn[0];  
  //  y = Inn[1];  
  double la = Inn[2]; // age
  double lv = Inn[3]; // 
  double IntL_val = Inn[4];
  double SuscepPopRatio = Inn[5];
  double VaccinePopRatio = Inn[6];
  double NormalisedVaccineRatio = Inn[7];  
  double VacFact = 0.05 + 0.95*VaccinePopRatio;

  double k5 = 2.3937;

  // cout << "VacFact: " <<  VacFact <<endl;
  // if(t<20)
  // { sd = 0.1 + t*0.2/20; } // f3 =0.89
  // else if(t<75)
  // { sd = 0.3; }  
  // else if(t<150)
  // { sd = 0.3 + (t-75.)*0.2/75; } // f3 =0.5  
  // else 
  // { sd = 0.5; } // f3 =0.5


       // plot [0:0.5] 1. -  1./(1.+ exp(-(x - 0.5)/0.1)) //for gnuplot

  //  f1 = 1./(1.+ exp(-(sigma - sigma_c)/bsigma));
  // f1 = 1.;
  // f2 = 1. -  1./(1.+ exp(-(health - health_c)/bhealth));
  // f2 = 1.;

  // double f5, f4 = a4*exp(-(la -b4)* (la -b4)/c4_sqr); // used till 23 May 2021
  double f5, f4;
  
  //Jun 2021, https://ncdc.gov.in/dashboard.php
  // if(la<=10./125.) {
  //   f4 = 0.0336; }
  // else if(la<=20./125.) {
  //   f4 = 0.0841; }
  // else if(la<=30./125.){
  //   f4 = 0.218; }
  // else if(la<=40./125.){
  //   f4 = 0.2193; }    
  // else if(la<=50./125.){
  //   f4 = 0.1722; }    
  // else if(la<=60./125.){
  //   f4 = 0.1405; }    
  // else if(la<=70./125.){
  //   f4 = 0.0866; }    
  // else if(la<=80./125.){
  //   f4 = 0.0354; }    
  // else if(la<=90./125.){
  //   f4 = 0.0092; }    
  // else {
  //   f4 = 0.0012;}
  

 if (TDatabase::ParamDB->PBE_P10<0.)
  {
  if(la<=11./125.) {
    f4 = 2.56919365672*TDatabase::ParamDB->PBE_P5; 
    }
  else if(la<=17./125.) {
    f4 = 4.70419796961*TDatabase::ParamDB->PBE_P6; 
    }
  else if(la<=44./125.){
    f4 =  1.04537706409*TDatabase::ParamDB->PBE_P7;
     }
  else if(la<=59./125.){
    f4 = 1.88168560379*TDatabase::ParamDB->PBE_P8;
     }    
  else {
    f4 = 0.42765485306*TDatabase::ParamDB->PBE_P9;   
    }  
  f4 *=NormalisedVaccineRatio;      
  // cout << "NormalisedVaccineRatio: " << NormalisedVaccineRatio << endl;   
  }   
 else
 {
    f4 = NormalisedVaccineRatio/5.0; // N_AgeGroup = 5;
 }

  // 
  double lvsim2 = (lv-vsymp)*(lv-vsymp);
  if(lv<vsymp) {
    f5 = k5*exp( -lvsim2/(2.*(vsymp/3.)*(vsymp/3.)) ); }
   else  {
    f5 = k5*exp( -lvsim2/(VacFact*2*((1.-vsymp)/3.)*((1.-vsymp)/3.)) );  }
  
  double norm = TDatabase::ParamDB->FS_U;

  B_Nuc =SuscepPopRatio*R_0*f1*sd*f4*f5*IntL_val/norm; 

   if(B_Nuc<0)
   {
    cout << "B_Nuc: " << B_Nuc << " IntL_val "  <<IntL_val << " : " << sd << 
     " : " << f5 << endl;
    // exit(0);
   }  

  Out[0] = 1; // growth is 1
  Out[1] = B_Nuc; 
  Out[2] = 0; //L_Max Bound value,  if DIRICHLET 
    // cout << "B_Nuc: " <<  Out[1] <<endl;
}


// Function to update the Internal Size parameters
 void GetLExampleData(int *L_NVertices, double *L_StartEnd, BoundCond1D **BoundConLminLMax, 
                      DoubleFunctND **GrowthAndB_Nuc, double **XPos)
  {
   // Internal Co-ordinates
   
   // Start of domain
   L_StartEnd[0] = 0.1;
   
   // End of Domain
   L_StartEnd[1] = 1.0;
   
   // Number of cells in domain
   L_NVertices[0] = 4;
   
   // Growth and Nucleation Kernels
   GrowthAndB_Nuc[0] = GrowthAndB_Nuc_L0;
   
   // Boundary Conditions
   BoundConLminLMax[0] = BoundCondition_LminLMax_L0;
   
   // Predefined X-Arrays, If you have any fixed x-arrays, those can be provided here
   XPos[0] = nullptr;

   }




// Function to obtain internal Nodal points
void GetInternalNodalPts(TFESpace1D *FeSpace_Intl, int &N_Intl_Levels, double *&IntlPosL)
{
  int i, j, k, l, m, r, N_Cells, N_RootPts, *RootPtIndex;
  int N_DOFs, N_LocalDOFs, N_Points;
  int N_AllLocalPoints;

  double L, L0, *xi, *eta, *L_loc, *L_loc_origOrder;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;

  Coll = FeSpace_Intl->GetCollection();
  N_Cells = Coll->GetN_Cells();

  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FeSpace_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    N_AllLocalPoints +=N_Points;
  }

  L_loc = new double [N_AllLocalPoints];
  L_loc_origOrder = new double [N_AllLocalPoints];
  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FeSpace_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);

    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, X, Y, AbsDetjk);

    for(j=0; j<N_Points; j++)
    {
      L_loc[N_AllLocalPoints] = X[j];
      N_AllLocalPoints++;
    }
  } // for(i=0; i<N_Cells; i++)

  memcpy(L_loc_origOrder, L_loc,  N_AllLocalPoints*SizeOfDouble);

  for(i=0; i<N_AllLocalPoints-1; i++)
   for(j=i; j<N_AllLocalPoints; j++)
    if(L_loc[i]> L_loc[j])
     {
      L= L_loc[i];
      L_loc[i]= L_loc[j];
      L_loc[j]= L;
     }

  L  = L_loc[0];
  N_RootPts = 1;

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      N_RootPts++;
      L = L_loc[i];
     }
   }

  IntlPosL= new double[N_AllLocalPoints];
  IntlPosL[0] = L_loc[0];
  N_RootPts = 1;
  L  = L_loc[0];

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      IntlPosL[N_RootPts] = L_loc[i];
      N_RootPts++;
      L = L_loc[i];
     }
   }

  delete [] L_loc;

  RootPtIndex = new int[N_AllLocalPoints];

  // find the index for the local points in the root points
  for(i=0; i<N_AllLocalPoints; i++)
   {
    L = L_loc_origOrder[i];
    l=0;
    r=N_RootPts;

    m = N_RootPts/2;
    L0 = IntlPosL[m];

    while(fabs(L-L0) > 1.e-8 )
     {
      if(L < L0)  //poin lies in left side
       {
        r = m;
       }
      else
       {
        l=m;
       }

      m= (l+r)/2;
      L0 = IntlPosL[m];
     } //  while ( 

    RootPtIndex[i] = m;
   }

  FeSpace_Intl->SetIntlPtIndexOfPts(RootPtIndex);
  FeSpace_Intl->SetN_RootNodalPts(N_RootPts);

  N_Intl_Levels = N_RootPts;

//  cout << N_AllLocalPoints << " N_RootPts  "  << N_RootPts << endl;
//   for(i=0; i<N_RootPts; i++)
//   cout << i << " L: "  << IntlPosL[i] << endl;
// exit(0);
}

// Thivin - Function to obtain the Quadrature points of the physical space
void GetPhysicalSpaceQuadPts(int &N_PhySpacePts, TFESpace3D *FESpace3D, double *&IntlX, double *&IntlY, double *&IntlZ)
{
  int i,j,k,l, m;
  int N_Cells, PolynomialDegree, N_Points;

  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D], AbsDetjk[MaxN_QuadPoints_3D];
  double *xi, *eta, *zeta, *weights;

  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId;
  TFE3D *Element;
  BF3DRefElements RefElement;
  QuadFormula3D QuadFormula;
  TQuadFormula3D *qf3;
  TRefTrans3D *F_K;

  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();

   m = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace3D->GetFE3D(i, cell);
    Element = TFEDatabase3D::GetFE3D(FEId);
    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
    PolynomialDegree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);

     switch(RefElement)
    {
      case BFUnitHexahedron:
        QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)F_K)->SetCell(cell);
        break;

      case BFUnitTetrahedron:
        QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)F_K)->SetCell(cell);
        break;
    }                                             // endswitch

//     cout << "QuadFormula: " << QuadFormula << endl;
    qf3 = TFEDatabase3D::GetQuadFormula3D(QuadFormula);
    qf3->GetFormulaData(N_Points, weights, xi, eta, zeta);

    if(i==0)
    {
      N_PhySpacePts = N_Points*N_Cells;

      IntlX = new double[N_PhySpacePts];
      IntlY = new double[N_PhySpacePts];
      IntlZ = new double[N_PhySpacePts];      
    }

    switch(RefElement)
    {
      case BFUnitHexahedron:
        ((THexaAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, AbsDetjk);
        break;

      case BFUnitTetrahedron:
        ((TTetraAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, AbsDetjk);
        break;
    }     

    for(j=0; j<N_Points; j++)
    {
      IntlX[m] = X[j];
      IntlY[m] = Y[j];
      IntlZ[m] = Z[j];
      m++;
    }
       
  } //   for(i=0; i<N_Cells; i++)
  
//      cout<< N_PhySpacePts << " N_PhySpacePts " << m <<endl;
//     for(i=0; i<N_PhySpacePts; i++)
//      cout<< i << " IntlX " << IntlX[i] << " IntlY " << IntlY[i]  << " IntlZ " << IntlZ[i]  <<  endl;
//    exit(0);
//   cout << " GetPhysicalSpaceQuadPts " << endl;
//   exit(0);
} // GetPhysicalSpaceQuadPts
