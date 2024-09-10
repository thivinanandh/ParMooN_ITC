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


// 1D Mesh Generation 

void Generate1DMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, z;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  
  N_V = N_Cells+1; 

//     double *X;
//   X = new double[N_V];
// 
//   for(i=0; i<N_V; i++)
//    X[i] = 1. + (1. - Start)*(tanh(2.75*(double(i)/double(N_Cells) - 1.)))/tanh(2.75);
// 
//   h = (End - Start)/N_Cells;
//   for(i=1; i<N_V; i++)
//    X[i] =  h*(double)i;
  
  
  if(N_Cells!=5)
   {
    printf("Error: N_Cells_Internal != 112, Generate1DMesh : %d \n", N_Cells); 
#ifdef _MPI
     MPI_Finalize();
#endif       
    // exit(0);
   } 

  
//   double X[113] = {0.001250000000000 , 0.004948299443979 , 0.006217671606001 , 0.007823152131089 , 
// 0.009849848192817 , 0.012405802335839 , 0.015627666201516 , 0.019687946143196 , 0.024804199643075 , 
// 0.028393316830700 , 0.031250666632101 , 0.035772828811153 , 0.039372952752786 , 0.042413137150944 , 
// 0.045070619540390 , 0.047446991851855 , 0.049606547407142 , 0.053436976294678 , 0.056785220390605 , 
// 0.059779281606086 , 0.062500166625532 , 0.065002773507971 , 0.067326227648337 , 0.069499512571455 , 
// 0.071544767303250 , 0.075317060474946 , 0.078745170560444 , 0.081898265580090 , 0.084825640943857 , 
// 0.087563813936174 , 0.090140678210433 , 0.092578078038834 , 0.094893477609017 , 0.097101084797404 , 
// 0.099212631825087 , 0.103185287607544 , 0.106873553596882 , 0.110323438897215 , 0.113570087453620 , 
// 0.116641021427289 , 0.119558244391400 , 0.122339657205393 , 0.125000041585273 , 0.130005277376165 , 
// 0.134652203946564 , 0.138998789266748 , 0.143089312023620 , 0.146958448032182 , 0.150633920105047 , 
// 0.154138281947656 , 0.157490157382677 , 0.163796361297538 , 0.169651123546876 , 0.175127479279515 , 
// 0.180281216202297 , 0.185156023145274 , 0.189786828693602 , 0.194202048758084 , 0.198425147902268 , 
// 0.206370468208258 , 0.213747007445228 , 0.220646784186760 , 0.227140086575032 , 0.233281959112378 , 
// 0.239116409077369 , 0.244679238288390 , 0.250000010253906 , 0.260010487342228 , 0.269304345055468 , 
// 0.277997519564343 , 0.286178568401423 , 0.293916843310072 , 0.301267789998808 , 0.308276515941197 , 
// 0.314980268830742 , 0.327592680129367 , 0.339302207508500 , 0.350254921410786 , 0.360562397349922 , 
// 0.370312013057423 , 0.379573625756073 , 0.388404067306966 , 0.396850266867541 , 0.412740909664792 , 
// 0.427493989953308 , 0.441293544971590 , 0.454280151067003 , 0.466563897289198 , 0.478232798228372 , 
// 0.489358457546173 , 0.500000002278646 , 0.520020957831925 , 0.538608674401517 , 0.555995024386395 , 
// 0.572357122891389 , 0.587833673431567 , 0.602535567444792 , 0.616553019893863 , 0.629960526177828 , 
// 0.655185349642304 , 0.678604405120686 , 0.700509833534509 , 0.721124785936175 , 0.740624017806563 , 
// 0.759147243604363 , 0.776808127061630 , 0.793700526500832 , 0.825481812641652 , 0.854987973672328 , 
// 0.882587084092698 , 0.908560296613240 , 0.956465591475152 , 1.000000000000000}; 

     //neu grid von ARAM, ref Bulk_3d4d.C
    //  double X[94] = {0, 0.000488281, 0.000615194, 0.000775097, 0.000976561, 0.00123039, 0.00155019, 0.00195312,0.00246078, 
    //                       0.00310039, 0.00390625, 0.00447154, 0.00492156, 0.0053016, 0.00563379, 0.00593084, 0.00620079,
    //                       0.00667959, 0.00709813, 0.00747238, 0.0078125, 0.00841576, 0.00894308, 0.00941462, 0.00984313, 
    //                       0.0106032, 0.0112676, 0.0118617, 0.0124016, 0.0133592, 0.0141962, 0.0149448, 0.015625, 0.0168315, 
    //                       0.0178862, 0.0188292, 0.0196863, 0.0204745, 0.0212064, 0.0218909, 0.0225352, 0.0231445, 0.0237233, 
    //                       0.0242752, 0.0248031, 0.0267184, 0.0283925, 0.0298896, 0.03125, 0.0336631, 0.0357723, 0.0393725, 
    //                       0.0450703, 0.0496062, 0.0534368, 0.056785, 0.0597791, 0.0625, 0.0673261, 0.0715446, 0.0753169,
    //                       0.0787451, 0.0818982, 0.0848255, 0.0875637, 0.0901406, 0.0948934, 0.0992125, 0.106873, 0.11357,
    //                       0.119558, 0.125, 0.134652, 0.143089, 0.150634, 0.15749, 0.169651, 0.180281, 0.189787, 0.198425, 
    //                       0.213747, 0.22714, 0.239117, 0.25, 0.269304, 0.286179, 0.301268, 0.31498, 0.360562, 0.39685, 
    //                       0.5, 0.629961, 0.793701, 1};

    double X[5] = {0, 0.25, 0.5, 0.75, 1};


  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_V]; 

  y=0.;
  z=0.;
  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
    Vetrex[i] = new TVertex(X[i], y, z);
   }

  Vetrex[N_Cells] = new TVertex(X[N_V-1], y, z);

  CellTree = new TBaseCell*[N_Cells];

   for (i=0;i<N_Cells;i++)
   {
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);
    ((TMacroCell *) CellTree[i])->SetSubGridID(0);

//  cout<< " x " <<CellTree[i]->GetN_Edges()<<endl;;
//  cout<< " x " <<TDatabase::RefDescDB[S_Line]->GetN_OrigEdges()<<endl;;

   }


//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
//     exit(0);

   Domain->SetTreeInfo(CellTree, N_Cells);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   // start joint(vertex)
   Joint = new TJointEqN(CellTree[0]);
   CellTree[0]->SetJoint(0, Joint);


   for(i=1;i<N_Cells;i++)
    {
     Joint = new TJointEqN(CellTree[i-1], CellTree[i]);

     CellTree[i-1]->SetJoint(1, Joint);
     CellTree[i]->SetJoint(0, Joint);
   } // for(i=0;i<N_Cells;i++)

   // end joint(vertex)
   Joint = new TJointEqN(CellTree[N_Cells-1]);
   CellTree[N_Cells-1]->SetJoint(1, Joint);

  delete []  Lines;
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
