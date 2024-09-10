// =======================================================================
//
// Purpose:     main program for Radiative Transfer Equation equation with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 27.04.19
// =======================================================================
 
#include <Domain.h>
#include <Database.h>
#include <SystemPBE3D.h>
#include <FEDatabase2D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <Output3D.h>
#include <LinAlg.h>

#include <MainUtilities.h>
#include <TimeDiscRout.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>


#include <stdlib.h>
#include<MacroCell.h>

// Thivin
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <ADISystem1D.h>

#ifdef _SMPI
#include "mpi.h"
#include <MeshPartition.h>
#endif

double bound = 0;
double timeC = 0;
// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TCD_3D/rte.h"
// #include "../Examples/TCD_3D/rteSpatial.h"
#include "../Examples/TCD_3D/Siminhale_PBE.h"

int main(int argc, char* argv[])
{
// ======================================================================
// variable declaration
// ======================================================================
  int i, j, l, m, N_SubSteps, ORDER, LEVELS, mg_level, N_Cells, N_DOF, img=1;
  int N_Active, mg_type, N_DOF_Solid, N_Solid_Levels, N_DoF_Intl;
  
  double *sol, *oldsol, *sol_tmp, *rhs, *oldrhs, t1, t2, errors[5], **Sol_array, **Rhs_array;
  double *Sol_SpatialSurface, *oldsol_all;
  double tau, end_time, *defect, *olderror, *olderror1, *H1T, *L2T, hmin, hmax;
  double start_time, stop_time,	 start_vtk=0, end_vtk=0, total_vtk=0, l2, H1, L2error_Max=0., L2error_Max_t;
	double start_assembling=0, end_assembling=0, total_assembling=0;
	double start_solve=0, end_solve=0, total_solve=0;
	double start_int=0, end_int=0, total_int=0;
  double *InternalPts, *Mat_Intl, val, *InternalScaling;  

  TDomain *Domain, *Domain_Solid;
  TDatabase *Database = new TDatabase();

  bool UpdateStiffnessMat, UpdateOnlyRhs,  ConvectionFirstTime;
  
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO, *PRM, *GEO_SOLID, *PRM_SOLID;
  char Name[] = "name";
  char Description[] = "description";
  char CString[] = "C";
  double Linfty=0;
  char SubID[] = "";
  
  int profiling;
#ifdef _SMPI
  int rank, size, out_rank;
  int MaxCpV, MaxSubDomainPerDof;
  
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  double time1, time2, temp, reduced_errors[4];

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;

  int Refine;
  int metisType[2] = {0,0};
#endif
  
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TFEDatabase2D *FEDatabase_surface = new TFEDatabase2D();
  TCollection *coll;
  TFESpace3D **Scalar_FeSpaces, *fesp[1], **Scalar_Solid_FeSpaces;
  TFEFunction3D *Scalar_FeFunction, **Scalar_FeFunctions, *Scalar_Solid_FeFunction, **Scalar_Solid_FeFunctions;
  TOutput3D *Output;
  TSystemPBE3D *System_Space;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = {D000, D100, D010, D001};

  std::ostringstream os;
  os << " ";   
    
  mkdir(vtkdir, 0777);

// ======================================================================
// set the database values and generate mesh
// ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  

 #ifdef _SMPI
  if(rank==0)
  { TDatabase::ParamDB->Par_P0=1; }
  else
  { TDatabase::ParamDB->Par_P0=0;  }  
  #endif

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  #ifdef _SMPI
  if(TDatabase::ParamDB->Par_P0)
  #endif
   {
    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();
   }

  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh  
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  PsBaseName = TDatabase::ParamDB->BASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  /* meshgenerator */
  if(TDatabase::ParamDB->MESH_TYPE==0)
    {Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {Domain->GmshGen(TDatabase::ParamDB->GEOFILE); }//gmsh mesh
  // else if(TDatabase::ParamDB->MESH_TYPE==2)   
    // {
      // Domain->TetrameshGen(TDatabase::ParamDB->GEOFILE);
    // } //tetgen mesh
    else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     }
  

  LEVELS = TDatabase::ParamDB->LEVELS;

  if(TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
  {
    TDatabase::ParamDB->UNIFORM_STEPS += (LEVELS-1);
    LEVELS = 1;
  }
  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();

  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
  #ifdef _SMPI
  if(TDatabase::ParamDB->Par_P0)
  #endif    
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }

//=========================================================================
// set data for multigrid
//=========================================================================  
  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;

  if(TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
     mg_type=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SCALAR = mg_type;
    }
  
  if(mg_type)
   {
    mg_level =  LEVELS + 1;
    ORDER = -1;
   }
  else
   {
    mg_level = LEVELS;
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
   }
  
  if(TDatabase::ParamDB->SOLVER_TYPE==GMG)
   {
  #ifdef _SMPI
  if(TDatabase::ParamDB->Par_P0)
  #endif
    {
    OutPut("=======================================================" << endl);
    OutPut("======           GEOMETRY  LEVEL ");
    OutPut(LEVELS << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level << "              ======" << endl);
    OutPut("=======================================================" << endl);   
    }
   }

  Scalar_FeSpaces = new TFESpace3D*[mg_level];  
  Scalar_FeFunctions = new TFEFunction3D*[mg_level]; 
  Sol_array = new double*[mg_level];
  Rhs_array = new double*[mg_level];
  
//=========================================================================
// construct all finite element spaces
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================

  for(i=0;i<LEVELS;i++)
   {   
    if(i)
     { Domain->RegRefineAll(); }
     
     coll=Domain->GetCollection(It_Finest, 0);
      
     // fespaces for scalar equation     
    Scalar_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);     

     //multilevel multigrid disc
    if(i==LEVELS-1 && mg_type==1) 
    {
      ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
      Scalar_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);
    } //  if(i==LEVELS-1 && i!=mg_level-1) 

//======================================================================
// construct all finite element functions
//======================================================================
  N_DOF = Scalar_FeSpaces[i]->GetN_DegreesOfFreedom();
  sol = new double[N_DOF];
  rhs = new double[N_DOF];
  Sol_array[i] = sol;
  Rhs_array[i] = rhs;   
         
  if(i==LEVELS-1 && mg_type==1) 
    {  
    N_DOF = Scalar_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    Sol_array[mg_level-1] = sol;
    Rhs_array[mg_level-1] = rhs;

    Scalar_FeFunction = new TFEFunction3D(Scalar_FeSpaces[mg_level-1], CString, CString, sol, N_DOF);   
    Scalar_FeFunctions[mg_level-1] = Scalar_FeFunction;
    }//   if(i==LEVELS-1 && mg_type==1) 
  else
    {
      Scalar_FeFunction  = new TFEFunction3D(Scalar_FeSpaces[i], CString, CString, sol, N_DOF);  
      Scalar_FeFunctions[i] = Scalar_FeFunction;
    }

#ifdef _MPI
     N_Cells = coll->GetN_Cells();
     printf("rank=%d\t N_Cells   : %d\t Dof all   :%d\n",rank,N_Cells,N_DOF);
#endif
   }// for(i=0;i<LEVELS;i++)

   oldrhs = new double[N_DOF];
   oldsol = new double[N_DOF];
   
#ifndef _MPI   
   N_Cells = coll->GetN_Cells();
  #ifdef _SMPI
  if(TDatabase::ParamDB->Par_P0)
  #endif
  {   
   OutPut("[INFO] N_Cells   : " << N_Cells <<endl);
   OutPut("[INFO] Dof all   : " << N_DOF  << endl);  
   OutPut(endl);
  }
#endif


//=========================================================================
// Spatial position at which the internal system is solverd
//=========================================================================
 double *PosX;
 PosX = new double[3*N_DOF]; 
 InternalScaling = new double[N_DOF];
 double *OldIntlRhs = new double[N_DOF];
 double fact;
 TFEVectFunct3D *Pos_VectFunction = new TFEVectFunct3D(Scalar_FeSpaces[mg_level-1], "Pos", "Pos", 
                                    PosX, N_DOF, 3);  
 Pos_VectFunction->GridToData();

// --- THIVIN ---
int N_PhySpacePts = N_DOF;


// Generate a Velocity FESpace for obtaining the Velocity values of bent pipe
TFESpace3D*  velocity_space = new TFESpace3D(coll, "Velocity Space", "Velocity Space", BoundCondition_Velocity, ORDER);

// Setup memory for velocity array
int N_VeloDof = velocity_space->GetN_DegreesOfFreedom();
double* velocity = new double [3*N_VeloDof];                                  
memset(velocity, 0, 3*N_VeloDof*SizeOfDouble);

// Generate FE Function 3D
TFEVectFunct3D* u = new TFEVectFunct3D(velocity_space, "U", "U", velocity, N_VeloDof, 3);
TFEFunction3D* u1 = u->GetComponent(0);
TFEFunction3D* u2 = u->GetComponent(1);
TFEFunction3D* u3 = u->GetComponent(2);

// Read the Solution from the Binary File
std::string u1_file = "Solution_u_000003.bin";
std::string u2_file = "Solution_v_000003.bin";
std::string u3_file = "Solution_w_000003.bin";

std::ifstream u1_in(u1_file, std::ios::in | std::ios::binary);
std::ifstream u2_in(u2_file, std::ios::in | std::ios::binary);
std::ifstream u3_in(u3_file, std::ios::in | std::ios::binary);

if (!u1_in.is_open() || !u2_in.is_open() || !u3_in.is_open())
{
  std::cerr << "Error: Cannot open the file " << std::endl;
  exit(1);
}

// Read the Solution from the Binary File
u1_in.read((char*)velocity, sizeof(double)* N_VeloDof);
u2_in.read((char*)velocity + sizeof(double)* N_VeloDof, sizeof(double)* N_VeloDof);
u3_in.read((char*)velocity + 2*sizeof(double)* N_VeloDof, sizeof(double)* N_VeloDof);

u1_in.close();
u2_in.close();
u3_in.close();

cout << "[INFO] Velocity Files successfully read from files " << endl;


// Generate 1D mesh for internal co-ordinates
TDomain *Domain_Intl = new TDomain(argv[1]);  
TCollection *Coll_Intl;
int N_Cells_Intl;
double L0, L1;
double* IntlPosL;

// Spatial position for 1D Data for Internal Co-ordinates
int N = 5;
L0 = 0.1;
L1 = 1.0;

Generate1DMesh(Domain_Intl, L0, L1, N);

Coll_Intl = Domain_Intl->GetCollection(It_Finest, 0);
N_Cells_Intl= Coll_Intl->GetN_Cells();

// Generate FE Space for Internal Co-ordinates
TFESpace1D*  FeSpace_Intl = new TFESpace1D(Coll_Intl, "Internal-1D Space", "Internal-1D Space", TDatabase::ParamDB->TEST_ORDER);


int N_Intl_Dof =  FeSpace_Intl->GetN_DegreesOfFreedom();
int N_Intl_Levels = N_Intl_Dof;

cout << "[INFO] N_Intl_Dof " << N_Intl_Dof << endl;
cout << "[INFO] N_Intl_Levels " << N_Intl_Levels << endl;

// N_Intl_Levels will be changed accordingly to 
// the no of points which are needed to evaluate internal nodal functionals
GetInternalNodalPts(FeSpace_Intl, N_Intl_Levels, IntlPosL);

// Generate layer Mass
double *LayerMass = new double[N_Intl_Levels];
for (i=0; i<N_Intl_Levels; i++)
{
  LayerMass[i] = pow((IntlPosL[i]), 3.0);
}

// Setup 1D Structures
TSquareStructure1D* SqStruct1D_Intl = new TSquareStructure1D(FeSpace_Intl);
SqStruct1D_Intl->Sort();
TSquareMatrix1D* M_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
TSquareMatrix1D* A_Intl = new TSquareMatrix1D(SqStruct1D_Intl);

cout << "[INFO] SqStruct1D_Intl->GetN_Rows() " << SqStruct1D_Intl->GetN_Rows() << endl;

// Obtain the Physical Space points
double* IntlX, *IntlY, *IntlZ;
double* XNodalPts =  new double[3 * N_DOF];
TFEVectFunct3D* PhySpaceNodalPts = new TFEVectFunct3D(Scalar_FeSpaces[mg_level-1], "PhySpaceNodalPts", "PhySpaceNodalPts", 
                                              XNodalPts, N_DOF, 3);

// Obtain the Physical Space Points
PhySpaceNodalPts->GridToData();
// transfer them to a double pointer array
IntlX = XNodalPts;
IntlY = XNodalPts + N_DOF;
IntlZ = XNodalPts + 2*N_DOF;

cout <<"[INFO] Physical Space Points obtained " << endl;


// Setup a System matrix to solve for internal Co-ordinates at all Physical Points in the actual 3D domain
TADISystem1D** ADI_System = new TADISystem1D*[N_DOF];

// Setup common solution variables for the ADI System
double* Sol_IntlLoc = new double[N_Intl_Dof]();
double* OldSol_IntlLoc = new double[N_Intl_Dof]();
double* B_IntlLoc = new double[N_Intl_Dof]();
double* defect_IntlLoc = new double[N_Intl_Dof]();

// Setup FeFunction for Internal Co-ordinates
TFEFunction1D* Intl_FeFunction = new TFEFunction1D(FeSpace_Intl, "Internal Co-ordinates", "Internal Co-ordinates", 
                                                Sol_IntlLoc, N_Intl_Dof);

// Setup the ADI System
for (int physical_quad = 0; physical_quad < N_DOF; physical_quad++)
{                          
  ADI_System[physical_quad] = new TADISystem1D(Intl_FeFunction, M_Intl, A_Intl, 3, XNodalPts,  Sol_IntlLoc, OldSol_IntlLoc, B_IntlLoc, defect_IntlLoc, NULL);
}

// System Matrices - Setup
cout << "[INFO] System Matrices Setup for Internal Co-ordinates" << endl;

// Setup solution arrays for internal points
double* Sol_Intl = new double[N_DOF * N_Intl_Dof]();
double* OldSol_Intl = new double[N_DOF * N_Intl_Dof]();
double* Rhs_Intl = new double[N_DOF * N_Intl_Dof]();
double* IntlSol_All_SUM = new double[N_Intl_Dof]();
double* IntlSol_All_SUM_T = new double[N_Intl_Dof]();
double* IntlSol_tmp2 = new double[N_Intl_Dof]();


// Set up memory for the internal system
int N_U = N_DOF;
int PhyDofTimesIntlLevels = N_Intl_Levels*N_U;
  
double* SolPbe = new double[PhyDofTimesIntlLevels];
double* B_Pbe = new double[PhyDofTimesIntlLevels];
  

Output = new TOutput3D(2, N_Intl_Levels, 1, 1, Domain);
Output->AddFEVectFunct(u);
Output->AddFEFunction(Scalar_FeFunction);
os.seekp(std::ios::beg);
Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str()); 

//======================================================================
// System_Space construction and solution
//======================================================================  

    /** interpolate the initial value */
    // Scalar_FeFunction->Interpolate(InitialCondition);   
    
    /** Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
     *  Solver: AMG_SOLVE (or) GMG  (or) DIRECT */   
    System_Space = new TSystemPBE3D(mg_level, Scalar_FeSpaces, Sol_array, Rhs_array,
                                  TDatabase::ParamDB->DISCTYPE, TDatabase::ParamDB->SOLVER_TYPE, 
                                  Scalar_FeFunctions, GetKernel, GetRhs, Domain_Intl);  

    
 #ifdef _SMPI
     if(rank==0)
 #endif
    printf("System_Space constructed\n");
  //      MPI_Finalize();
  //  exit(0);
    /** initilize the system  in space with the functions defined in Example file */
    // last argument aux is used to pass  additional fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass it with NULL 
    System_Space->Init(BilinearCoeffs, BoundCondition, BoundValue, BoundCondition, BoundValue, NULL,
                       SurfBoundCondition, SurfBoundValue);
    exit(0);
  //Get internal points
  InternalPts=System_Space->GetIntlPoints();
  N_DoF_Intl= System_Space->GetN_IntlPoints();
  Sol_SpatialSurface = new double[N_DOF*N_DoF_Intl];
  oldsol_all = new double[N_DOF*N_DoF_Intl];
  memset(Sol_SpatialSurface, 0, N_DOF*N_DoF_Intl*SizeOfDouble);     


 //======================================================================
 // produce outout at t=0
 //======================================================================
  // VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
  // Output = new TOutput3D(2, 2, 1, 1, Domain);
  // Output->AddFEFunction(Scalar_FeFunction);

  // if(TDatabase::ParamDB->WRITE_VTK)
  //  {
  //   os.seekp(std::ios::beg);
  //   if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
  //   else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
  //   else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
  //   else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
  //   else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
  //   Output->WriteVtk(os.str().c_str());
  //   img++;
  //   }   

  // interpolate the initial condition 
  System_Space->Interpolate(Exact, Sol_SpatialSurface); 
       

  /** measure errors to known solution */
  double oldl2_RTE, l2L2_RTE = 0.;
  if(TDatabase::ParamDB->MEASURE_ERRORS)
   {  
    System_Space->GetErrors(Exact_RTE, Scalar_FeFunction, Sol_SpatialSurface, l2);

  #ifdef _SMPI
  if(TDatabase::ParamDB->Par_P0)
  #endif
    {
    if(L2error_Max<l2)
      {
       L2error_Max=l2;
       L2error_Max_t=TDatabase::TimeDB->CURRENTTIME;
      } 

    OutPut("Time: " << TDatabase::TimeDB->CURRENTTIME);
    OutPut(" L2: " << l2);
    OutPut(" L_Infty " << L2error_Max << endl);
    oldl2_RTE = l2;
    }

   } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  


   
  coll->GetHminHmax(&hmin,&hmax);
 
#ifdef _MPI
  MPI_Allreduce(&hmax, &temp, 1, MPI_DOUBLE, MPI_MAX, Comm);
  hmax = temp;
  temp = 0;
#endif

#ifdef _SMPI
   if(rank == out_rank)
#endif
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

 TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;
//======================================================================
// time disc loop
//======================================================================    
// parameters for time stepping scheme
  m = 0;
  N_SubSteps = GetN_SubSteps();
  end_time = TDatabase::TimeDB->ENDTIME; 
   
 // assemble the internal system matrix
 // timestep is neede to assemble, so set time parameters
  SetTimeDiscParameters(1);
  // System_Space->AssembleIntlMat(); 
  Mat_Intl = System_Space->GetIntlMat();

  if (TDatabase::TimeDB->TIME_DISC>2)
   { 
    cout<< "Only Euler and CN schemes are implemented for RTE" <<endl;
    exit(0);
   }

  // Assemble M, A and Rhs are t=0
  System_Space->AssembleABRhs(); 
  System_Space->AssembleRhs();  // to copy rhs to old rhs
  UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
  UpdateOnlyRhs = TRUE; //check BilinearCoeffs in example file

  /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
  memcpy(oldrhs, rhs, N_DOF*SizeOfDouble);    
  

// =================================================================================================
/** time loop starts */
// =================================================================================================
while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {
   m++;
   TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

   for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
    {
     SetTimeDiscParameters(1);

#ifdef _SMPI
     if(rank == out_rank)
#endif
     if(m==1)
      {
       OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
       OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
       OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
       OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
      }
 
     tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     TDatabase::TimeDB->CURRENTTIME += tau;

#ifdef _SMPI
  if(rank == out_rank)
#endif
     OutPut(endl << "CURRENT TIME: "<<TDatabase::TimeDB->CURRENTTIME << endl);  

  if(UpdateStiffnessMat || UpdateOnlyRhs)
   {  
    if(UpdateOnlyRhs)
     { System_Space->AssembleRhs(); } // since f deepends on t
    else
     { System_Space->AssembleABRhs();  }
   }
   else
   {
    UpdateOnlyRhs = TRUE; 
   }

  // Strang splitting L-update t^{n+1/2}
   for(i=0;i<N_DoF_Intl;i++)
   {
     sol_tmp = Sol_SpatialSurface+(i*N_DOF);
    //  TDatabase::TimeDB->CURRENTTIME -= tau/2.;
     System_Space->AssembleIntlMat(i, tau, PosX, fact, InternalScaling);       //solve internal system
     Dsum(N_DOF, 1, fact, InternalScaling, sol_tmp, sol_tmp); /** z := alpha*x + beta*y */
    //  TDatabase::TimeDB->CURRENTTIME += tau/2.;
   }

   /** copy the sol(t^{n+1/2}) as old sol */  
   memcpy(oldsol_all, Sol_SpatialSurface, N_DOF*N_DoF_Intl*SizeOfDouble);  

    // Strang splitting X-update t^{n+1}
   for(i=0;i<N_DoF_Intl;i++)
   {
     sol_tmp = Sol_SpatialSurface+(i*N_DOF);
     oldsol = oldsol_all+(i*N_DOF);
     /**  M:= M + (tau*TDatabase::TimeDB->THETA1)*A */
     /** rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs + [ M - (tau*THETA2)A]*oldsol **/
    //  System_Space->AssembleSystMat(oldrhs, oldsol, rhs, sol_tmp,0,0,0     
	   System_Space->AssembleSystMat(oldrhs, oldsol, rhs, sol_tmp, 
                    InternalPts[3*i],  InternalPts[(3*i)+1],  InternalPts[(3*i)+2], SystemRhs
 #ifdef _MPI
                                       , Rhs_array
 #endif 
                                      );
     // Solve the system
     System_Space->Solve(sol_tmp);

     //restore mass matrix
     System_Space->RestoreMassMat();
   } // for(i=0;i<N_DoF_Int

  //  // Strang splitting L-update t^{n+1}
  //  for(i=0;i<N_DoF_Intl;i++)
  //  {
  //    sol_tmp = Sol_SpatialSurface+(i*N_DOF);
  //    System_Space->AssembleIntlMat(i, tau/2., PosX, fact, InternalScaling);       //solve internal system
  //    Dsum(N_DOF, 1, fact, InternalScaling, sol_tmp, sol_tmp); /** z := alpha*x + beta*y */
  //  }
  } // for(l=0;l<N_SubSteps;l++) 
// //======================================================================
// // produce output
// //======================================================================
// #ifdef _SMPI
//     if(profiling){
//       total_vtk += (end_vtk-start_vtk);
//       start_vtk = MPI_Wtime();
//     }
//   if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
//     if(TDatabase::ParamDB->WRITE_VTK)
//       Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
//         img++;

//     if(profiling)	end_vtk = MPI_Wtime();
// #else
//     if(profiling){
//       total_vtk += (end_vtk-start_vtk);
//       start_vtk = GetTime();
//     }
   
//    if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
//     if(TDatabase::ParamDB->WRITE_VTK)
//      {
//       os.seekp(std::ios::beg);
//        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//             else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//       Output->WriteVtk(os.str().c_str());
//       img++;
//      }   
//     if(profiling)	end_vtk = GetTime();
// #endif


// cout<<"tau " << TDatabase::TimeDB->TIMESTEPLENGTH << " : " << tau << endl;
//======================================================================
// measure errors to known solution
//======================================================================    
   if(TDatabase::ParamDB->MEASURE_ERRORS)
    {    

    // for(i=0; i<N_DoF_Intl; ++i)
    //  {
    //   TDatabase::ParamDB->P11 = InternalPts[3*i];
    //   TDatabase::ParamDB->P12 = InternalPts[3*i+1];
    //   TDatabase::ParamDB->P13 = InternalPts[3*i+2];
    //   Scalar_FeFunction->Interpolate(Exact);
      
    //   memcpy(Sol_SpatialSurface+(i*N_DOF), sol, N_DOF*SizeOfDouble); 
    // } 

     System_Space->GetErrors(Exact_RTE, Scalar_FeFunction, Sol_SpatialSurface, l2);

  #ifdef _SMPI
  if(TDatabase::ParamDB->Par_P0)
  #endif
    {
     if(L2error_Max<l2)
      {
       L2error_Max=l2;
       L2error_Max_t=TDatabase::TimeDB->CURRENTTIME;
      } 

     OutPut("Time: " << TDatabase::TimeDB->CURRENTTIME);
     OutPut(" L2: " << l2);
     OutPut(" L_Infty " << L2error_Max << endl);
 
      if(m>1)
        {
          l2L2_RTE += (l2*l2 + oldl2_RTE *oldl2_RTE)*0.5*TDatabase::TimeDB->TIMESTEPLENGTH;
          
          OutPut("L2 : " << l2 <<   " L2_Max: " << L2error_Max << " L2_Max_Time: " <<  L2error_Max_t <<" CurrentTime: ");
          OutPut(TDatabase::TimeDB->CURRENTTIME <<  " l2(0,T;L2) " << sqrt(l2L2_RTE) << endl);
         }

      oldl2_RTE = l2;
    }
 
   } //  if(TDatabase::ParamDB->MEASURE_ERRORS) 
// exit(0);

 } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

  CloseFiles();
 
#ifdef _SMPI
 MPI_Finalize();
#endif


} // end main
