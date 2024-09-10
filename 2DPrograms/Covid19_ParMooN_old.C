// =======================================================================
//
// Purpose:     main program for Covid-19 model in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 04.04.2020

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemTCD2D.h>
#include <SystemTCD2D_ALE.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>

#include <ADISystem.h>
#include <ADISystem1D.h>
#include <LineAffin.h>
#include <NodalFunctional1D.h>
#include <SquareMatrix1D.h>
#include <SquareStructure1D.h>
#include <BaseCell.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TCD_2D/exp.h"
// #include "../Examples/TCD_2D/SinCos1.h"
// #include "../Examples_All/TCD_2D/Time3.h"
// #include "../Examples/TCD_2D/exp_0.h"
//    #include "../Examples/TCD_2D/exp_2.h"
// #include "../Examples_All/TCD_2D/exp_1.h"
#include "../MainUsers/Sashi/TCD_2D/covid19-example.h"

// =======================================================================
// include file for internal functions
// =======================================================================
#include "Covid19.h"

int main(int argc, char* argv[])
{
  int i, j, l, m, N_SubSteps, N_Cells, N_XDOF, img=1, N_G;
  int N_Active;

  double *sol, *rhs, *oldrhs, t1, t2, errors[5], Linfty;
  double tau, end_time, *defect, olderror, olderror1, hmin, hmax;
  
  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
  char *VtkBaseName;
  const char vtkdir[] = "VTK"; 
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll;
  TFESpace2D *Scalar_FeSpace, *fesp[1];
  TFEFunction2D *Scalar_FeFunction;
  TOutput2D *Output;
  TSystemTCD2D *SystemMatrix;  
  TAuxParam2D *aux;
  MultiIndex2D AllDerivatives[3] = {D00, D10, D01}; 
  
  std::ostringstream os;
  os << " ";   
  
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  // set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based 
  Domain = new TDomain(argv[1]);  
  TDomain *Domain_Intl = new TDomain();
  TDomain *Domain_Intl_l2 = new TDomain();

  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 2;
   OpenFiles();

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  
  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */ 
  if(TDatabase::ParamDB->MESH_TYPE==0)
    {
      Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); 
      OutPut("PRM-GEO used for meshing !!!" << endl);
    } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {
       Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
       OutPut("GMSH used for meshing !!!" << endl);
    }//gmsh mesh
  else if(TDatabase::ParamDB->MESH_TYPE==2)    //triangle mesh
     {  
      OutPut("Triangle.h used for meshing !!!" << endl);
      TriaReMeshGen(Domain);
     } 
    else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     }
      
  // refine grid up to the coarsest level
  for(i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
#ifdef __LOCKDOWNMODEL__  
   Domain->RefineallxDirection();
#else
    Domain->RegRefineAll();  
#endif

  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain->PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(vtkdir, 0777);

  //=========================================================================
  // order of finite element spaces
  //=========================================================================
#ifdef __LOCKDOWNMODEL__  
TDatabase::ParamDB->ANSATZ_ORDER = 0;
#endif

  int PBE_ORDER = TDatabase::ParamDB->PBE_P8;
  int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  // int VELOCITYORDER = TDatabase::ParamDB->VELOCITY_SPACE;

  // initializ time
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH;
  SetTimeDiscParameters();

  double t3 = GetTime();
  double total_time = t3 - total_time;
  // SetPolynomialDegree();

  // check the example file, to activate
  int VeloFunction=TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD;

  BoundCond cond_Lmin, cond_Lmax;
  BoundCondition_LminLMax(cond_Lmin, cond_Lmax);

  //=====================================================================================
  // construct fespace  for population balance equation in physical space
  // FeSpace and memory for all matrices
  // assume that the convection and the reaction coefficient terms (if any) are independent
  // of internal coordinates, so lhs matrices are same for all lelvels of internal coordinate
  //=====================================================================================
  //=====================================================================================
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells (space) : " << N_Cells <<endl);
  
//   TBaseCell *cell;
//   double x, y, x_mean;
//   int N_Joints;
//   for(i=0;i<N_Cells;i++)
//    {
//     x_mean = 0.;
//     cell=coll->GetCell(i);
//     N_Joints = cell->GetN_Edges();

//     for(j=0;j<N_Joints;j++)
//      { 
//       (cell->GetVertex(j))->GetCoords(x, y);
//       x_mean +=x;
//      }
//       cout << "x_mean " << x_mean/(double)N_Joints << endl;
//    }
// exit(0);

  // fespaces for scalar equation 
  Scalar_FeSpace = new TFESpace2D(coll, (char*)"fe space", (char*)"solution space", 
                                          BoundCondition_Psd, ORDER, NULL);
   
  N_XDOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
  N_Active =  Scalar_FeSpace->GetActiveBound();
  OutPut("DOF PBE spatial space : "<< N_XDOF  << endl);
  OutPut("DOF active            : "  << N_Active << endl);

//======================================================================
// construct all finite element functions
//======================================================================
  sol = new double[N_XDOF];
  rhs = new double[N_XDOF];
  oldrhs = new double[N_XDOF];
    
  memset(sol, 0, N_XDOF*SizeOfDouble);
  memset(rhs, 0, N_XDOF*SizeOfDouble);

  Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, (char*)"sol", (char*)"sol", sol, N_XDOF); 
  
  //interpolate the initial value
  Scalar_FeFunction->Interpolate(InitialCondition_Psd);

#ifndef __LOCKDOWNMODEL__
//======================================================================
// SystemMatrix construction and solution
//======================================================================    
  // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
  // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
  SystemMatrix = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);
  
  // initilize the system matrix with the functions defined in Example file
  SystemMatrix->Init(BilinearCoeffs_Psd, BoundCondition_Psd, BoundValue_Psd);
     
  // assemble the system matrix with given aux, sol and rhs 
  // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
  // otherwise, just pass with NULL 
  SystemMatrix->AssembleMRhs(NULL, sol, rhs);   
#endif

//=========================================================================
// Construct FeSpace and all data for the internal domain
// FESpace is same in internal coordinate for all  QuadPts
//=========================================================================
  double *IntlX, *IntlY;
  int MaxN_PtsForNodal, N_Xpos=0;

 #ifdef __OSNODALPT__
    GetNodalPtsADI_System(N_Xpos, Scalar_FeSpace, MaxN_PtsForNodal, IntlX, IntlY); 
 #else
     GetQuadPtsADI_System(N_Xpos, Scalar_FeSpace, IntlX, IntlY); 
 #endif

  OutPut("Dof PBE Internal space N_Xpos : "<<  N_Xpos << endl);

  // interlanl l_1 mesh
  int N = int(TDatabase::ParamDB->REACTOR_P11);
  double L0 = TDatabase::ParamDB->REACTOR_P12;
  double L1 = TDatabase::ParamDB->REACTOR_P13;

  Generate1DMesh(Domain_Intl, L0, L1, N);
  TCollection *Coll_Intl = Domain_Intl->GetCollection(It_Finest, 0);
  int N_Cells_Intl= Coll_Intl->GetN_Cells();

  cout<< "N_Cells_Internal " << N_Cells_Intl <<endl;

  // Finite difference in internal 
  if(TDatabase::ParamDB->INTL_DISCTYPE==FD)
   { 
    TDatabase::ParamDB->TEST_ORDER = 1;
   }

  char IString[] = "I";
  TFESpace1D *FeSpace_Intl = new TFESpace1D(Coll_Intl, IString, IString, TDatabase::ParamDB->TEST_ORDER);
  int dGDisc, N_LDof, N_L1nodal;
  N_LDof = FeSpace_Intl->GetN_DegreesOfFreedom();

  if(TDatabase::ParamDB->TEST_ORDER<-9)
   {
    FeSpace_Intl->SetAsDGSpace(); 
    dGDisc = 1;
    TDatabase::ParamDB->INTL_DISCTYPE=DG;  // no supg method
   }

  double *IntlPosL;
  // nodal points to evaluate nodal values
  GetInternalNodalPts(FeSpace_Intl, N_L1nodal, IntlPosL);
  OutPut("Dof PBE Internal space : "<<  N_L1nodal << endl);
  TSquareStructure1D *SqStruct1D_Intl = new TSquareStructure1D(FeSpace_Intl);
  SqStruct1D_Intl->Sort();
  TSquareMatrix1D *S_Intl, *M_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
  TSquareMatrix1D *K_Intl, *A_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
  if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
   {
    S_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
    K_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
   }

 double *Sol_IntlLoc = new double[N_L1nodal];
 double *OldSol_IntlLoc = new double[N_L1nodal];
 double *B_IntlLoc = new double[N_L1nodal];
 double *defect_IntlLoc = new double[N_L1nodal];
 double *Sol_Xpos_IntOmegaL2L1 = new double[N_Xpos];

 memset(Sol_IntlLoc, 0, N_L1nodal*SizeOfDouble);
 memset(OldSol_IntlLoc, 0, N_L1nodal*SizeOfDouble);
 memset(B_IntlLoc, 0, N_L1nodal*SizeOfDouble);
 memset(defect_IntlLoc, 0, N_L1nodal*SizeOfDouble);

 TFEFunction1D *FeFunction_Intl = new TFEFunction1D(FeSpace_Intl, IString, IString, Sol_IntlLoc,  N_L1nodal);
//=========================================================================
// l_2 internal domain
//=========================================================================
  N = int(TDatabase::ParamDB->REACTOR_P15);

  Generate1DMesh(Domain_Intl_l2, L0, L1, N);
  TCollection *Coll_Intl_l2 = Domain_Intl_l2->GetCollection(It_Finest, 0);
  int N_Cells_Intl_l2= Coll_Intl_l2->GetN_Cells();

  cout<< "N_Cells_Internal_l2 " << N_Cells_Intl_l2 <<endl;
  TFESpace1D *FeSpace_Intl_l2 = new TFESpace1D(Coll_Intl_l2, IString, IString, TDatabase::ParamDB->TEST_ORDER);
  int N_L2Dof = FeSpace_Intl->GetN_DegreesOfFreedom();

  if(TDatabase::ParamDB->TEST_ORDER<-9)
   {
    FeSpace_Intl_l2->SetAsDGSpace();    
   }

  double *IntlPosL_l2;
  int N_L2nodal;
  GetInternalNodalPts(FeSpace_Intl_l2, N_L2nodal, IntlPosL_l2);
  OutPut("Dof PBE Internal_l2 space : "<< N_L2nodal << endl);
  TSquareStructure1D *SqStruct1D_Intl_l2 = new TSquareStructure1D(FeSpace_Intl_l2);
  SqStruct1D_Intl_l2->Sort();
  TSquareMatrix1D *S_Intl_l2, *M_Intl_l2 = new TSquareMatrix1D(SqStruct1D_Intl_l2);
  TSquareMatrix1D *K_Intl_l2, *A_Intl_l2 = new TSquareMatrix1D(SqStruct1D_Intl_l2);
  if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
   {
    S_Intl_l2 = new TSquareMatrix1D(SqStruct1D_Intl_l2);
    K_Intl_l2 = new TSquareMatrix1D(SqStruct1D_Intl_l2);
   }

//  TADISystem1D **ADI_System_l2 = new TADISystem1D*[N_Xpos];
 double *Sol_IntlLoc_l2 = new double[N_L2nodal];
 double *OldSol_IntlLoc_l2 = new double[N_L2nodal];
 double *B_IntlLoc_l2 = new double[N_L2nodal];
 double *defect_IntlLoc_l2 = new double[N_L2nodal];

 memset(Sol_IntlLoc_l2, 0, N_L2nodal*SizeOfDouble);
 memset(OldSol_IntlLoc_l2, 0, N_L2nodal*SizeOfDouble);
 memset(B_IntlLoc_l2, 0, N_L2nodal*SizeOfDouble);
 memset(defect_IntlLoc_l2, 0, N_L2nodal*SizeOfDouble);

 TFEFunction1D *FeFunction_Intl_l2 = new TFEFunction1D(FeSpace_Intl_l2, IString, IString, Sol_IntlLoc_l2,  N_L2nodal);

 int N_l1ADI_Systems = N_Xpos*N_L2nodal;
 int N_l2ADI_Systems = N_Xpos*N_L1nodal;
 //=========================================================================
 // l_1 system
 //=========================================================================

 TADISystem1D **ADI_System = new TADISystem1D*[N_l1ADI_Systems];

 cout << " ADI_System[i] loop is wrong, (ii is missing), see Covid19_L3_ParMooN code" <<endl;
  exit(0);

 for(int ii=0; ii<N_L2nodal; ii++)
  for(i=0; i<N_Xpos; i++)
  {
   ADI_System[i] = new TADISystem1D(FeFunction_Intl, M_Intl, A_Intl, S_Intl, K_Intl, IntlX[i], IntlY[i],
                                    InitialCondition_Psd_Intl,  BoundValue_LMin, BoundValue_LMax,
                                    Sol_IntlLoc, OldSol_IntlLoc, B_IntlLoc, defect_IntlLoc, IntlPosL,
#ifdef __L1GROWTHNUCL__   
                                    Get_GrowthAndB_Nuc 
#else
                                     NULL
#endif
                                     , IntlPosL_l2[ii]
                                     );
   }
     
//=========================================================================
// l_2 system
//=========================================================================
 TADISystem1D **ADI_System_l2 = new TADISystem1D*[N_l2ADI_Systems];

 for(int ii=0; ii<N_L1nodal; ii++)
  for(i=0; i<N_Xpos; i++)
  {
   ADI_System_l2[i] = new TADISystem1D(FeFunction_Intl_l2, M_Intl_l2, A_Intl_l2, S_Intl_l2, K_Intl_l2, IntlX[i], IntlY[i],
                                    InitialCondition_Psd_Intl,  BoundValue_LMin_l2, BoundValue_LMax_l2,
                                    Sol_IntlLoc_l2, OldSol_IntlLoc_l2, B_IntlLoc_l2, defect_IntlLoc_l2, IntlPosL_l2,
#ifdef __L2GROWTHNUCL__   
                                     Get_GrowthAndB_Nuc_L2 
#else
                                     NULL
#endif
                                     , IntlPosL[ii]
                                     );
   }

//=========================================================================
// memory allocate for all vectors and construction of PBE fefunction in
// physical space (for all internal levels)
//=========================================================================
int N_NodalDofPos_All = N_Xpos*N_L1nodal*N_L2nodal;

double *IntL1L2value = new double[N_Xpos];
double *Sol_L2nodalXposL1nodal = new double[N_NodalDofPos_All];
double *OldSol_L1nodalXposL1nodal = new double[N_NodalDofPos_All];
double *PBE_IntlPtValues = new double[N_NodalDofPos_All];
double *PBE_IntlPtRhsValues = new double[N_NodalDofPos_All];
// bool *DirichletBDPt = new bool[N_Xpos*N_L2nodal];
// bool *DirichletBDPt_l2 = new bool[N_Xpos*N_L1nodal];

int N_Dof_All = N_Xpos*N_LDof*N_L2nodal;
double *Sol_L2nodalXposL1dof =  new double[N_Dof_All];
double *Sol_XposL1dof, *Sol_L1dof;
double *Rhs_L2nodalXposL1dof = new double[N_Dof_All];
double *Rhs_XposL1dof, *Rhs_L1dof;
double *OldPBE_IntlPtRhsValuesT = new double[N_Dof_All];
// double *Oldrhs_IntlLoc = new double[N_LDof;
// double *Oldrhs_IntlLoc_l2 = new double[N_L2Dof];

memset(OldPBE_IntlPtRhsValuesT, 0, N_Dof_All*SizeOfDouble);
memset(Rhs_L2nodalXposL1dof, 0, N_Dof_All*SizeOfDouble);

N_Dof_All = N_L2nodal*N_L1nodal*N_XDOF;
double *Sol_L2nodalL1nodalXdof = new double[N_Dof_All];
memset(Sol_L2nodalL1nodalXdof, 0, N_Dof_All*SizeOfDouble);    
double *B_Pbe = new double[N_Dof_All];
memset(B_Pbe, 0, N_Dof_All*SizeOfDouble);    
double *RhsArray_Pbe = new double[N_Dof_All];
memset(RhsArray_Pbe, 0, N_Dof_All*SizeOfDouble);
double *OldRhsArray_Pbe = new double[N_Dof_All];
memset(OldRhsArray_Pbe, 0, N_Dof_All*SizeOfDouble);

// set output parameters
int start_pbe, end_pbe, Out_Level = N_L1nodal/4;

  if(cond_Lmin==DIRICHLET && !dGDisc)
   { start_pbe = 1; }
  else
   { start_pbe = 0; }

  if(cond_Lmax==DIRICHLET && !dGDisc) { end_pbe = N_L1nodal-1; }
  else { end_pbe = N_L1nodal; }

 double *SolPbe_Output = new double[N_XDOF];
 memset(SolPbe_Output, 0, N_XDOF*SizeOfDouble);
 TFEFunction2D *PbeFunctions = new  TFEFunction2D(Scalar_FeSpace, (char*)"I", (char*)"I", SolPbe_Output, N_XDOF);

// prepare output
VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    

Output = new TOutput2D(0,  1, 1, 1, Domain);
Output->AddFEFunction(PbeFunctions);

os.seekp(std::ios::beg);
Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());

int N_l1l2Dofs = N_L2nodal*N_Xpos;
TFEFunction2D **ScalarPbeFunctions = new TFEFunction2D*[N_l1l2Dofs];
   
TFEVectFunct2D *Pbe_xdir = new TFEVectFunct2D(Scalar_FeSpace, (char*)"I", (char*)"I", Sol_L2nodalL1nodalXdof, N_XDOF, N_l1l2Dofs);
  
for(j=0; j<N_l1l2Dofs; j++)
   ScalarPbeFunctions[j] = Pbe_xdir->GetComponent(j);

// interpolate initial solution at L direction nodal-interpolation point levels
TDatabase::TimeDB->CURRENTTIME=0.;

memset(Sol_L2nodalXposL1nodal, 0, N_NodalDofPos_All*SizeOfDouble);
double *Sol_XposL1nodal,  *Sol_L1nodalXdof, *Sol_Xdof;

for(int ii=0; ii<N_L2nodal; ii++)
{
 Sol_XposL1nodal = Sol_L2nodalXposL1nodal+(ii*N_L1nodal*N_Xpos);

 for(i=0;i<N_Xpos;i++)
  {
   Sol_L1nodalXdof = Sol_XposL1nodal + i*N_L1nodal;
  
   ADI_System[i]->Interpolate(Sol_L1nodalXdof, InitialCondition_Psd_Intl);
  }
} //for(int ii=0; ii<N_

// exit(0);


/** transfer sol from internal to spatial */
//  ScalarPbeFunctions_low= ScalarPbeFunctions if PBE_ORDER == ORDER
for(int ii=0; ii<N_L2nodal; ii++)
{
 Sol_XposL1nodal = Sol_L2nodalXposL1nodal+(ii*N_L1nodal*N_Xpos);
 Sol_L1nodalXdof = Sol_L2nodalL1nodalXdof+(ii*N_L1nodal*N_XDOF);

 #ifdef __OSNODALPT__
 GetSolFromNodalPtVales(N_L1nodal, MaxN_PtsForNodal, ScalarPbeFunctions, Sol_L1nodalXdof,
                        N_XDOF, Sol_XposL1nodal, N_Xpos);
 #else
 GetSolFromQuadPtVales(N_L1nodal, MaxN_PtsForNodal, ScalarPbeFunctions, Sol_L1nodalXdof, 
                       N_XDOF, Sol_XposL1nodal, N_Xpos);
 #endif

    //  for(j=0; j<N_XDOF; j++)
    // { 
    //  for(i=0; i<N_L1nodal; i++)
    //   cout<< i << " Sol_L1nodalXdof " << Sol_L1nodalXdof[j*N_L1nodal + i] <<endl;
    // cout<<endl; 
    // }
    // cout<<endl; 

}

 double l2 = 0.;
 double H1 = 0.;
 if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
   for(int ii=0; ii<N_L2nodal; ii++)
    {
     Sol_L1nodalXdof = Sol_L2nodalL1nodalXdof+(ii*N_L1nodal*N_XDOF);
     TDatabase::ParamDB->REACTOR_P30 = IntlPosL_l2[ii];
  
     GetOSError(N_L1nodal, ScalarPbeFunctions,  FeSpace_Intl, Sol_L1nodalXdof, N_XDOF, errors, Exact_Psd_Intl);
     l2 += errors[0];
     H1 += errors[1];
    }

   l2 /=(double)N_L2nodal;
   H1 /=(double)N_L2nodal;

   //small l2
   errors[5] = (TDatabase::TimeDB->TIMESTEPLENGTH)*l2;

   OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
   OutPut(" L2: " << sqrt(l2));
   OutPut(" H1-semi: " << sqrt(H1) << endl);
  }



if(TDatabase::ParamDB->WRITE_VTK)
for(j=0; j<N_L1nodal; j++)
 {
  TDatabase::ParamDB->REACTOR_P29=IntlPosL[j];
  memcpy(SolPbe_Output, Sol_L2nodalL1nodalXdof+j*N_XDOF, N_XDOF*SizeOfDouble);

  os.seekp(std::ios::beg);
  if(img<10) os <<"VTK/"<<  VtkBaseName<<j<<".0000"<<img<<".vtk" << ends;
  else if(img<100) os <<"VTK/"<<  VtkBaseName<<j<<".000"<<img<<".vtk" << ends;
  else if(img<1000) os <<"VTK/"<<  VtkBaseName<<j<<".00"<<img<<".vtk" << ends;
  else if(img<10000) os <<"VTK/"<<  VtkBaseName<<j<<".0"<<img<<".vtk" << ends;
  else  os <<"VTK/"<<  VtkBaseName<<j<<"."<<img<<".vtk" << ends;
  Output->WriteVtk(os.str().c_str());
 } 


//======================================================================
// prepare at t=0
//======================================================================
 
// assemble the mass matrix for the internal space once for all time steps
  if(TDatabase::ParamDB->INTL_DISCTYPE==FD)
   { AssembleFD_M_StartPtDirichlet(M_Intl); }
  else
   { AssembleM_StartPtDirichlet(FeSpace_Intl, M_Intl);}
 
   /** store the M Mat values, so that no need to restore after every sptl pt. calculatio*/
   int N_Entries = M_Intl->GetN_Entries();
   double *MatValues = new double[N_Entries];
   memcpy(MatValues,  M_Intl->GetEntries(),  N_Entries*SizeOfDouble);
  
  // TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;

  double GrowthFact_X, GrowthParameters[4];
  double C_old = TDatabase::ParamDB->REACTOR_P2;
  double C_new = TDatabase::ParamDB->REACTOR_P2;
  //======================================================================
  // time disc loop
  //======================================================================    
  // parameters for time stepping scheme
   m = 0;
   N_SubSteps = GetN_SubSteps();
   end_time = TDatabase::TimeDB->ENDTIME; 

   UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
   UpdateRhs = FALSE; //check BilinearCoeffs in example file
   ConvectionFirstTime=FALSE;

   // time loop starts
   while(TDatabase::TimeDB->CURRENTTIME < end_time)
   {
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(l=0; l<N_SubSteps; l++) // sub steps of fractional step theta
    {
      SetTimeDiscParameters(1);

      if(m==1)
      {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
      }

      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
       

      GrowthFact_X = (1. - tau*TDatabase::TimeDB->THETA2*C_old)/(1. + tau*TDatabase::TimeDB->THETA1*C_new);

      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME << " GrowthFact_X: " << GrowthFact_X << endl);  

// OutPut(C_old << " GrowthFact_X: " << C_new << endl);  
// exit(0);

      for(int ii=0; ii<N_L2nodal; ii++)
       {
        Sol_L1nodalXdof = Sol_L2nodalL1nodalXdof+(ii*N_L1nodal*N_XDOF);

        #ifdef __LOCKDOWNMODEL__
         for(int i=0; i<N_L1nodal; i++)
          {
           Sol_Xdof = Sol_L1nodalXdof+i*N_XDOF;
           
           for(int jj=0; jj<N_XDOF; ++jj)
              {
               Sol_Xdof[jj] *=GrowthFact_X;
              }
          }
        #else

        cout <<"NOT __LOCKDOWNMODEL__ Not yet implemented " <<endl;
        exit(0);

          //copy rhs to oldrhs
         memcpy(oldrhs, rhs, N_XDOF*SizeOfDouble); 

        // unless the stiffness matrix or rhs change in time, it is enough to 
        // assemble only once at the begning
        if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
         {
          SystemMatrix->AssembleARhs(NULL, sol, rhs);

          // M:= M + (tau*THETA1)*A
          // rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs +[M-(tau*THETA2)A]*oldsol
          // note! sol contains only the previous time step value, so just pass 
          // sol for oldsol
          SystemMatrix->AssembleSystMat(oldrhs, sol, rhs, sol);
          ConvectionFirstTime = FALSE;
        }
     
        // solve the system matrix 
        SystemMatrix->Solve(sol, rhs);
    
        // restore the mass matrix for the next time step    
        // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
        if(UpdateStiffnessMat || UpdateRhs)
         {
          SystemMatrix->RestoreMassMat();
         } 
         #endif
      } //   for(int ii=0; ii 

     //===============================================================================
     //PBE system X-direction solution -- end
     //PBE system l1-direction solution -- start  
     //===============================================================================
     GetSol_L2_X_L1_and_Sol_X_IntOmegaL1L2(FeSpace_Intl, FeSpace_Intl_l2, Scalar_FeSpace, N_Xpos, N_L1nodal, N_L2nodal, Sol_L2nodalL1nodalXdof, 
                                           Sol_L2nodalXposL1nodal, Sol_L2nodalXposL1dof, Sol_Xpos_IntOmegaL2L1, IntL1L2value);

     for(int ii=0; ii<N_L2nodal; ii++)
      {

       Sol_XposL1dof = Sol_L2nodalXposL1dof+(ii*N_LDof*N_Xpos);
       Rhs_XposL1dof = Rhs_L2nodalXposL1dof+(ii*N_LDof*N_Xpos);

       for(i=0;i<N_Xpos;i++)
       {
         Rhs_L1dof = Rhs_XposL1dof+ (i*N_LDof);
         Sol_L1dof = Sol_XposL1dof + (i*N_LDof);

         GrowthParameters[0] = IntlPosL[i];  // l1
         GrowthParameters[1] = IntlPosL_l2[ii]; //l2
         GrowthParameters[2] = IntL1L2value[i]; //\int_\Omega_(l1xl2)

          ADI_System[i]->SolveAllLevels(MatValues, Sol_L1dof, Rhs_L1dof, BilinearCoeffs_Psd_Intl, GrowthParameters, tau, cond_Lmin, cond_Lmax); 
       }

        Sol_XposL1nodal = Sol_L2nodalXposL1nodal+(ii*N_L1nodal*N_Xpos);
        Sol_L1nodalXdof = Sol_L2nodalL1nodalXdof+(ii*N_L1nodal*N_XDOF);          
        L_Sol2Nodal(FeSpace_Intl, Sol_XposL1dof, N_L1nodal, Sol_XposL1nodal, N_Xpos);   

       // transpose the solution, 
       #ifdef __OSNODALPT__     
          GetSolFromNodalPtVales(N_L1nodal, MaxN_PtsForNodal, ScalarPbeFunctions, Sol_L1nodalXdof, N_XDOF, Sol_XposL1nodal, N_Xpos);
       #else
         GetSolFromQuadPtVales(N_L1nodal, MaxN_PtsForNodal, ScalarPbeFunctions, Sol_L1nodalXdof, N_XDOF, Sol_XposL1nodal, N_Xpos);
       #endif
      } // for(int ii=0; ii<N_PBE_Intl
     //===============================================================================
     //PBE system l1-direction solution -- end
     //PBE system l2-direction solution -- start  
     //===============================================================================

      //===============================================================================
      //PBE system l2-direction solution -- end 
      //===============================================================================

    } // for(l=0;l<N_SubSteps;l++) 
    // //======================================================================
    // // produce outout
    // //======================================================================
    if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
    {
      GetSol_L2_X_L1_and_Sol_X_IntOmegaL1L2(FeSpace_Intl, FeSpace_Intl_l2, Scalar_FeSpace, N_Xpos, N_L1nodal, N_L2nodal, Sol_L2nodalL1nodalXdof, 
                                            Sol_L2nodalXposL1nodal, Sol_L2nodalXposL1dof, Sol_Xpos_IntOmegaL2L1, IntL1L2value);

     for(i=1; i<=N_Xpos; i++)
      OutPut("Covid Population in the state " <<  i << " : " << IntL1L2value[N_Xpos - i] <<endl);

    if(TDatabase::ParamDB->WRITE_VTK)
    for(j=0; j<N_L1nodal; j++)
    {
      TDatabase::ParamDB->REACTOR_P29=IntlPosL[j];
      memcpy(SolPbe_Output, Sol_L2nodalL1nodalXdof+j*N_XDOF, N_XDOF*SizeOfDouble);

      os.seekp(std::ios::beg);
    if(img<10) os <<"VTK/"<<  VtkBaseName<<j<<".0000"<<img<<".vtk" << ends;
      else if(img<100) os <<"VTK/"<<  VtkBaseName<<j<<".000"<<img<<".vtk" << ends;
      else if(img<1000) os <<"VTK/"<<  VtkBaseName<<j<<".00"<<img<<".vtk" << ends;
      else if(img<10000) os <<"VTK/"<<  VtkBaseName<<j<<".0"<<img<<".vtk" << ends;
      else  os <<"VTK/"<<  VtkBaseName<<j<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
 } 

    }
    //======================================================================
    // measure errors to known solution
    //======================================================================    
    // if(TDatabase::ParamDB->MEASURE_ERRORS)
    // {
    //   Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs_Psd, aux, 1, fesp, errors);


    //   OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
    //   OutPut(" L2: " << errors[0]);
    //   OutPut(" H1-semi: " << errors[1] << endl);

    //   errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
    //   olderror = errors[0];
    //   OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

    //   errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
    //   OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
    //   olderror1 = errors[1];   
      
      
    //   if(Linfty<errors[0])
    //   Linfty=errors[0];

    //   OutPut( "Linfty " << Linfty << endl);      
    // } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  

  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)
//======================================================================
// produce final outout
//======================================================================
  
// for(int ii=0; ii<N_L2nodal; ii++)
//    {
//     for(j=0; j<N_L1nodal ; j++)
//     { 
//      for(i=0; i<N_XDOF; i++)
//       cout<< i << " Sol_L1nodalXdof " << Sol_L2nodalL1nodalXdof[ii*N_L1nodal*N_XDOF + j*N_XDOF + i] <<endl;
//     cout<<endl; 
//     // exit(0);
//     }
//     cout<<endl; 
//     }

// // exit(0);


    //  if(TDatabase::ParamDB->WRITE_VTK)
    //   {
    //    os.seekp(std::ios::beg);
    //     if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
    //      else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
    //       else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
    //        else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
    //         else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
    //    Output->WriteVtk(os.str().c_str());
    //    img++;
    //   }
      
      
  CloseFiles();
  
  
  return 0;
} // end main









