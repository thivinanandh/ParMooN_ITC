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
#include <SystemRTE.h>
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
#include "../Examples/TCD_3D/rteSpatialInternal.h"

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
  TSystemRTE *System_Space;
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
   OutPut("N_Cells   : " << N_Cells <<endl);
   OutPut("Dof all   : " << N_DOF  << endl);  
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

//======================================================================
// System_Space construction and solution
//======================================================================  

    /** interpolate the initial value */
    // Scalar_FeFunction->Interpolate(InitialCondition);   
    
    /** Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
     *  Solver: AMG_SOLVE (or) GMG  (or) DIRECT */   
    System_Space = new TSystemRTE(mg_level, Scalar_FeSpaces, Sol_array, Rhs_array,
                                  TDatabase::ParamDB->DISCTYPE, TDatabase::ParamDB->SOLVER_TYPE, 
                                  Scalar_FeFunctions, GetKernel, GetRhs);  
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
