// =======================================================================
// Purpose:     main program for Covid model with new kernels of ParMooN
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 04.04.2020
// =======================================================================
#ifdef _MPI
# include "mpi.h"
#endif

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

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <SystemADI1D.h>
#include <SystemADI1D_3L.h>
// =======================================================================
// include current example
// =======================================================================
// #include "../MainUsers/Sashi/TCD_2D/covid19-example.h"
// #include "../MainUsers/Sashi/TCD_2D/covid19-1AugInit.h"
// #include "../MainUsers/Sashi/TCD_2D/covid19-Mar21.h"
#include "../MainUsers/Sashi/TCD_2D/covid19-Mar21-KA.h"
// #include "../MainUsers/Sashi/TCD_2D/covid19-Mar21-JH.h"
// #include "../MainUsers/Sashi/TCD_2D/covid19-Mar21-PY.h"
// #include "../MainUsers/Sashi/TCD_2D/covid19-Mar21-TN.h"

// =======================================================================
// include file for internal functions
// =======================================================================
#include "PopulationBalance.h"

int main(int argc, char* argv[])
{
  #ifdef _MPI
  const int root = 0;
  int rank, mpi_size;
  double t_par1, t_par2;
  char  name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);  
  #endif 

  int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF, img=1, N_G;
  int N_Active;

  double TotalPopulation, *sol, *rhs, *oldrhs, t1, t2, errors[5], Linfty;
  double tau, end_time, *defect, olderror, olderror1, hmin, hmax;
  double AntibodyPop, AntibodyPopTot, VaccinePop, VaccinePopTot;

  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
  char *VtkBaseName, *CovidPopulationName, *CovidNucleationName, *CovidRecoverdName, *CovidDeathName;
  const char vtkdir[] = "VTK"; 
  const char Populationdir[] = "PopulationData"; 

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
  
  int N_AgeGroup=5;
  double MaxAgeOfGroup[5];
  MaxAgeOfGroup[0] = 11./125.; 
  MaxAgeOfGroup[1] = 17./125.; 
  MaxAgeOfGroup[2] = 44./125.; 
  MaxAgeOfGroup[3] = 59./125.; 
  MaxAgeOfGroup[4] = 125./125.;
  int N_EffAntibdyDays = 350;
 
  // India
  // double StatePopulation_All [32] = {417036, 53903393,  1570458, 35607039, 124799926, 1158473, 29436231, 18710922, 
  //           1586250, 63872399, 28204692, 7451955, 13606320, 38593948, 67562686, 35699443, 
  //            289023, 50073183, 85358965, 123144223, 3091545, 46356334, 1413542, 30141373, 
  //            81032689, 690251, 77841267, 39362732, 4169794, 237882725,  11250858, 99609303};
 // LD changed from 73183 to 50073183

  // Karnataka - Distric-wise population
  double StatePopulation_All [32] = {2141784., 3096512., 5354714., 1139694., 13626081., 1903546.,
                                     1076085., 1364310., 1136558., 1803209., 2335183., 2102415., 
                                     2099469., 1157112., 1840223., 1759087., 2976434., 560797., 
                                     1695959., 1594147., 1851525., 3374961., 2187486., 1140574., 
                                     1864370., 2784099., 1307243., 1516250., 2571844., 1412445., 1., 1};
  //area sq.mi                                     
  double State_Area_ALL[32] = {2539, 1642, 5180, 872, 850, 2103, 1970,  1747, 2780, 3260, 1760, 1720, 1640, 1798,
                           2631, 1862, 4228, 1584, 1532, 2776, 1915, 2646, 3260, 1373, 3273, 4092, 1500,
                           3973, 2179, 2021, 10000, 10000};
  //CIR Wave-I
  double State_CIR[32] = {31., 49., 95., 46., 23., 64., 72., 28., 54., 85., 34., 62., 13., 14.,
                          44., 72., 60., 56., 59., 38., 54., 34., 76., 76., 31., 82., 34., 33.,
                          112., 62., 1., 1.};
  // //CIR Wave-I
  // double State_CIR2[32] = {, 1., 1.};
  // double PopulationDensity[32];
  int N_Districts = 30;
    // 270: 28 Mar 2020 from 1 Jul
  // 380: 15 JUl 2021
  int NewVariant_day = 270;  // 270: 28 Mar 2020 from 1 Jul
  int NewVariant_day_Wave3 = -1; // no new variant  
  // int NewVariant_day_Wave3 = 380;   // 380: 15 JUl 2021
  // int NewVariant_day_Wave3 = 472;   // 472: 15 Oct 2021  



  // //Jarkhand - Distric-wise population
  // double StatePopulation_All [32] = {2351056., 1188890., 1700963., 3060315., 1506444., 2615068.,
  //                                    1507974., 2787840., 1497448., 1168743., 1977324., 901788.,
  //                                    606349.,  816535.,  828755.,  526441.,  1026481., 2211451.,
  //                                    1082365., 3322248., 1311646., 1214164., 683519.,  1712665.,
  //                                    1., 1., 1., 1., 1., 1., 1., 1 };

  //  TN- Distric-wise population
  // double StatePopulation_All [40] ={875677, 2965243, 5390209, 4011332, 3022860, 1747938, 2505339, 
  //                                   2612023, 1589526, 1353025, 2169634, 1234812, 2180578, 3524372, 
  //                                   1873893, 2002857, 853057, 655659, 1877280, 1569996, 1403921, 
  //                                   4039185, 1553357, 1632847, 2790832, 1445243, 4324601, 1466561, 
  //                                   2030204, 3157856, 1931693, 1289702, 2875700, 2859255, 1872521, 
  //                                   2427883, 2253054, 1., 1., 1 };


  //PY - Distric-wise population
  // TDatabase::ParamDB->UNIFORM_STEPS = 1;
  // double StatePopulation_All[4] = {264093, 55155, 1253431, 73371};

  std::ostringstream os;
  os << " ";   
  
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  // set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based 
  Domain = new TDomain(argv[1]);  
  // #ifdef _MPI  
  // if(rank==0)
  // { TDatabase::ParamDB->Par_P0=1; }
  // else
  // { TDatabase::ParamDB->Par_P0=0;  }  
  // #endif

  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 2;
   OpenFiles();

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
  mkdir(vtkdir, 0777);
  mkdir(Populationdir, 0777);

  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */ 
  if(TDatabase::ParamDB->MESH_TYPE==0)
    {
      Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); 
#ifdef _MPI
 if(TDatabase::ParamDB->Par_P0==1)
 #endif  
      OutPut("PRM-GEO used for meshing !!!" << endl);
      
    } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {
       Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
       OutPut("GMSH used for meshing !!!" << endl);
    }//gmsh mesh
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
  // if(TDatabase::ParamDB->WRITE_PS)
  //   Domain->PS("Domain.ps", It_Finest, 0);

  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  double *StatePopulation;
  double *State_Area;
  int N_AllCell;
 #ifdef _MPI
  
  TotalPopulation = 0;
  for(i=0; i<N_Districts; i++)
    TotalPopulation +=StatePopulation_All[i];

  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;++i)
   coll->GetCell(i)->SetGlobalCellNo(i);

  N_AllCell = N_Cells;

  // cout << "N_AllCell : " << N_AllCell <<endl;
  double **CovidPopulation_tmp;
  double *CovidPopulation1, *CovidPopulation2, *CovidPopulation3;
  CovidPopulation_tmp = new double*[N_AgeGroup];

  if(rank==0)
  {
   for(i=0; i<N_AgeGroup; i++)
    CovidPopulation_tmp[i] = new double[N_AllCell];

  //  CovidPopulation_tmp1 = new double[N_Cells];
  //  CovidPopulation_tmp2 = new double[N_Cells];      
  }

  CovidPopulation1 = new double[N_AllCell];
  CovidPopulation2 = new double[N_AllCell];
  CovidPopulation3 = new double[N_AllCell];

 #ifdef __LOCKDOWNMODEL__ 
 if(N_Cells<mpi_size)
  {
   if(TDatabase::ParamDB->Par_P0==1)
    printf("N_Cells less than N_MPI_Processes: N_Cells = %d\n", N_Cells);
   MPI_Finalize();
   exit(0); 
  }
  
  int N_OwnCells = N_Cells/mpi_size; 
  int N_RemainingCells = N_Cells%mpi_size; 
 
  // if(TDatabase::ParamDB->Par_P0==1)
  // printf("N_OwnCells = %d, N_RemainingCells = %d\n", N_OwnCells, N_RemainingCells);
  
  int AdditionCell=0;
  if(rank<N_RemainingCells)
   AdditionCell = 1 ;

  N_Cells = N_OwnCells + AdditionCell;

  TBaseCell **OwnCells;
  OwnCells = new TBaseCell*[N_Cells];
  int *GLOB_cellIndex = new int[N_Cells];

  int *N_Cells_ALL, *DispArray;
  N_Cells_ALL = new int[mpi_size];

  MPI_Allgather(&N_Cells, 1, MPI_INT, N_Cells_ALL, 1, MPI_INT, MPI_COMM_WORLD);

  int disp = 0;
  for(i=0;i<rank;++i)
   disp +=N_Cells_ALL[i];

  for(i=0;i<N_Cells;++i)
  {
   OwnCells[i] = coll->GetCell(disp+i);
   GLOB_cellIndex[i] = disp+i;
  }

  // printf("Rank %d, disp = %d, N_Cells = %d\n", rank, disp, N_Cells);

  Domain->ReplaceTreeInfo(N_Cells, OwnCells, GLOB_cellIndex, N_Cells);
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

  // if(rank==0)
  // {
   DispArray = new int[mpi_size];

  //  for(i=0;i<mpi_size;++i)
  //  {
  //   DispArray[i] = 0;
  //  }
   DispArray[0] = 0;
  // printf("Rank %d, N_Cells = %d\n", rank, N_Cells);
   for(i=1;i<mpi_size;++i)
    {
     DispArray[i] = DispArray[i-1] + N_Cells_ALL[i-1];
    //  if(rank==0)
    //  printf("Rank %d, N_Cells = %d, Disp = %d \n", i, N_Cells_ALL[i-1],  DispArray[i]);
    }
  // }

  //starting index
  StatePopulation = StatePopulation_All+DispArray[rank];
  State_Area = State_Area_ALL+DispArray[rank];

  // printf("Rank %d, StatePopulation = %f \n", rank,  StatePopulation[0]);
  //  MPI_Finalize();
  //  exit(0);
 #endif 
 #endif 
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  // printf("Rank %d, N_Cells = %d\n", rank, N_Cells);

  // fespaces for scalar equation 
  Scalar_FeSpace = new TFESpace2D(coll, (char*)"fe space", (char*)"solution space", 
                                          BoundCondition_Spatial, ORDER, NULL);
   
  N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
  N_Active =  Scalar_FeSpace->GetActiveBound();

 #ifdef _MPI
  printf("Rank, %d,  Dof all      : %d \n",rank, N_DOF );
  printf("Rank, %d,  Dof active   : %d \n", rank, N_Active);
#else
  OutPut("dof all      : "<< setw(10) << N_DOF  << endl);
  OutPut("dof active   : "<< setw(10) << N_Active << endl);
#endif
  //  MPI_Finalize();
  //  exit(0);
  //======================================================================
  // construct all finite element functions
  //======================================================================
  sol = new double[N_DOF];
  rhs = new double[N_DOF];
  oldrhs = new double[N_DOF];
    
  memset(sol, 0, N_DOF*SizeOfDouble);
  memset(rhs, 0, N_DOF*SizeOfDouble);

  Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, (char*)"sol", (char*)"sol", sol, N_DOF); 
 
  //======================================================================
  // SystemMatrix construction and solution
  //======================================================================    
  // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
  // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
  SystemMatrix = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);
    
  // initilize the system matrix with the functions defined in Example file
  SystemMatrix->Init(BilinearCoeffs_Spatial, BoundCondition_Spatial, BoundValue_Spatial);
      
  // assemble the system matrix with given aux, sol and rhs 
  // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
  // otherwise, just pass with NULL 
  SystemMatrix->AssembleMRhs(NULL, sol, rhs);   

  //=========================================================================
  // Spatial position at which the internal system is solverd
  //=========================================================================
  double *Xpos;
  int MaxN_PtsForNodal, N_Xpos=0;

 #ifdef __OSNODALPT__
    GetNodalPts(N_Xpos, Scalar_FeSpace, MaxN_PtsForNodal, Xpos); 
 #else
    cout << "Not Yet Implemented for Covid!" <<endl;
    exit(0);

    GetQuadPts(N_Xpos, Scalar_FeSpace, Xpos); 
 #endif
 
  //======================================================================
  // Insternal space construction and solution 
  //======================================================================   
  // double PopulDensFactor = TDatabase::ParamDB->;
  // for(i=0; i<N_Xpos; i++)
  //  {
  //  PopulationDensity[i] = StatePopulation[i]/(State_Area[i]*PopulDensFactor);  
  //  }

  //  if(rank==0)
  //  for(i=0; i<N_Districts; i++)
  //   OutPut("PopulationDensity :: "  << StatePopulation_All[i]/(State_Area_ALL[i]*PopulDensFactor)<<endl);

  //  printf("Rank, %d,  N_Xpos : %d \n",rank, N_Xpos );
  //  MPI_Finalize();
  //  exit(0);

  int ii, N_ADISystems = 3, N_LnodalPos[3], N_LDof[3], L_NVertices[3];
  double *LnodalPos[3], L_StartEnd[6], **XPos;
  DoubleFunctND  *GrowthAndB_Nuc[3];
  TSystemADI1D *SystemADI[3];
  TSystemADI1D_3L *ADISystem3L;
  BoundCond1D *BoundConLminLMax[3];

  XPos = new double *[3];
  XPos[0] = nullptr;
  XPos[1] = nullptr;
  XPos[2] = new double[N_AgeGroup+1];

  GetLExampleData(L_NVertices, L_StartEnd, BoundConLminLMax, GrowthAndB_Nuc, XPos);

  ADISystem3L = new TSystemADI1D_3L(SystemADI, L_NVertices, L_StartEnd, BoundConLminLMax, GrowthAndB_Nuc, Scalar_FeFunction, XPos);


// MPI_Finalize();
// exit(0);

  for(ii=0; ii<N_ADISystems; ii++)
   {
    N_LnodalPos[ii] = SystemADI[ii]->GetN_NodalPts();
    LnodalPos[ii] = SystemADI[ii]->GetNodalPt_Coord();
    N_LDof[ii] = SystemADI[ii]->GetN_Dof();
   }

   //ld
  i= SystemADI[0]->GetN_Cells();
  int N_LDistXSum = i*N_LnodalPos[2]*N_LnodalPos[1];

  // OutPut(endl << "N_LDistXSum  :  "<< N_LDistXSum << endl); 
 
  //solution arrays
  int N_Nodals_All = N_Xpos*N_LnodalPos[2]*N_LnodalPos[1]*N_LnodalPos[0];  
  int N_LNodals_All = N_LnodalPos[2]*N_LnodalPos[1]*N_LnodalPos[0];  
  double *Sol_XdofLdof = new double[N_Nodals_All];
  double *React_L0L1,  *React_L0,  *ReactionFactor = new double[2*N_Nodals_All];
  double *Sol_XDofLDof = new double[N_DOF*N_LDof[2]*N_LDof[1]*N_LDof[0]];
  double *LDistXSum = new double[N_LDistXSum];
  bool DistSolOutput =FALSE;
  double *sum_receiv = new double[N_LDistXSum];
  double *LRecoDistXSum = new double[N_LDistXSum];
  double *LRecoDistXSum_temp = new double[N_LNodals_All];
  memset(LRecoDistXSum, 0, N_LDistXSum*SizeOfDouble);  
  memset(LRecoDistXSum_temp, 0, N_LNodals_All*SizeOfDouble);  

  int N_LnodalXdof = N_DOF*N_LnodalPos[2]*N_LnodalPos[1]*N_LnodalPos[0];  
  double *Sol_LdofXdof = new double[N_LnodalXdof];
  double *Recovered_Xdof, *Recovery_LdofXdof = new double[N_LnodalXdof];

  // initiaize internal systems
  ADISystem3L->Init(N_Xpos, Xpos, N_LDof, N_LnodalPos, LnodalPos, Sol_XdofLdof, N_AgeGroup);

  //======================================================================
  // parameters for time stepping scheme
  //======================================================================  
  m = 0;
  N_SubSteps = GetN_SubSteps();
  end_time = TDatabase::TimeDB->ENDTIME;

#ifdef _MPI
 if(TDatabase::ParamDB->Par_P0==rank)
 #endif  
  {
   OutPut(endl << "Total Population:  "<< TotalPopulation/1.e6  <<" million  "<<  endl); 
   OutPut(endl << "CURRENT TIME:  "<< TDatabase::TimeDB->CURRENTTIME << endl); 
  }
  
#ifdef __LOCKDOWNMODEL__
  double InitialPopulation[50];
#ifdef _MPI
  // if(N_AllCell!=32)
  //   {

  //    if(TDatabase::ParamDB->Par_P0==1)
  //     OutPut("Lock down model must have 32 States/UT"<< N_AllCell <<endl);
  //     //  MPI_Finalize();
  //     // exit(0);  
  //   }
#endif   
  char *DataFile;
  DataFile = new char[50];
  strcpy(DataFile, "InitialData.dat");

  ReadData(N_AllCell, DataFile, InitialPopulation);

  //Interpolate 
  ADISystem3L->Interpolate(coll, InitialPopulation, InitialValues);
 #else 
  //Interpolate 
  ADISystem3L->Interpolate(InitialValues);
 #endif

  // XnodalLnodal->XnodalLdof->LdofXdof (note XnodalLdof is needed for IntL)
  ADISystem3L->XnodalLnodalToLdofXdof(MaxN_PtsForNodal, Sol_LdofXdof);

  //  MPI_Finalize();
  //  exit(0); 
  //======================================================================
  // produce outout at t=0
  //======================================================================
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;   
  CovidPopulationName = TDatabase::ParamDB->MAPFILE;   
  CovidNucleationName = TDatabase::ParamDB->PODFILE;   
  CovidRecoverdName = TDatabase::ParamDB->POD_FILENAME; 
  CovidDeathName = TDatabase::ParamDB->SNAP_FILENAME; 
  char *CovidNucleationNameAge[5] = {"Nucleation_11Yrs", "Nucleation_18Yrs", 
                                     "Nucleation_45Yrs", "Nucleation_60Yrs",
                                     "Nucleation_Above60Yrs"}; 

  Output = new TOutput2D(2, 2, 1, 1, Domain);

  Output->AddFEFunction(Scalar_FeFunction);

  int i_la, i_ld, i_lv;
  if(TDatabase::ParamDB->WRITE_VTK)
   {
    for(j=0;j<N_LNodals_All;j++)
     {
      // i_la =  (j/N_LnodalPos[0])/N_LnodalPos[1];
      // i_ld = (j/N_LnodalPos[0])%N_LnodalPos[1];
      // i_lv = j%N_LnodalPos[0];
      // cout << i_la << " : " << i_ld << " : " << i_lv <<endl;

      memcpy(sol, Sol_LdofXdof +j*N_DOF, N_DOF*SizeOfDouble);
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<j<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<j<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<j<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<j<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<j<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
     }   
    img++;
   }

  //======================================================================
  // Set variables at t=0
  //======================================================================
  int jj;
  double *Sol_Xdof, GrowthFact_X=1., RecoveryRatio=0.;

  //assemble M matirces of system L
  for(ii=0; ii<N_ADISystems; ++ii)
   {
    SystemADI[ii]->AssembleMassMat();
   }

  // ld in conservative form, to include B_nuc
  // if(rank==0)
  SystemADI[0]->AssembleAdvectMat(TRUE);
  SystemADI[1]->AssembleAdvectMat(TRUE);  
  SystemADI[2]->AssembleAdvectMat(TRUE);

  double *IntValue = new double[N_Xpos];
  double *IntValue_old = new double[N_Xpos];
  double *RecoveredValue = new double[N_Xpos];
  double *RecoveredValue_temp = new double[N_Xpos];     
  double *RemovedValue = new double[N_Xpos];  
  double *NucValue = new double[N_Xpos];

  memset(RemovedValue, 0, N_Xpos*SizeOfDouble);  
  memset(RecoveredValue, 0, N_Xpos*SizeOfDouble);  

  double PopRatio[5], InfectRatio[5];
  PopRatio[0] = TDatabase::ParamDB->PBE_P0;
  PopRatio[1] = TDatabase::ParamDB->PBE_P1;
  PopRatio[2] = TDatabase::ParamDB->PBE_P2;
  PopRatio[3] = TDatabase::ParamDB->PBE_P3;
  PopRatio[4] = TDatabase::ParamDB->PBE_P4;

  InfectRatio[0] = TDatabase::ParamDB->PBE_P5;  
  InfectRatio[1] = TDatabase::ParamDB->PBE_P6;
  InfectRatio[2] = TDatabase::ParamDB->PBE_P7;
  InfectRatio[3] = TDatabase::ParamDB->PBE_P8 ;
  InfectRatio[4] = TDatabase::ParamDB->PBE_P9;
  int AcclerateVacRate = TDatabase::ParamDB->KDP_c_p;

  //New variant induction day
  NewVariant_day_Wave3=TDatabase::ParamDB->KDP_w_sat_2;

  double *B_NucValue[5], *B_NucValue_tmp[5], *B_NucValue_old[5];
  double *ImmunPopulation[5], *ImmunPopulation_ALL, *B_NucValue_ALL, B_NucValue_Factor;
  double *VaccinatedPopulation[5];
  double *SuscepPopRatio[15], **AntibodyPopulation[5], **AntibodyPopulation_All[5];
  double TempTotal, TempTotal_all, TempAgeTotal_all, *TempSum = new double[N_Xpos], *TempAgeSum = new double[N_AgeGroup] ;

  B_NucValue_ALL = new double [N_AllCell];
  ImmunPopulation_ALL = new double [N_Xpos];

  N_EffAntibdyDays = int(end_time);
  jj = N_EffAntibdyDays+int(end_time); // antibody remains for 250 days
  for(i=0;i<N_AgeGroup; i++)
   {
    B_NucValue[i]= new double[N_Xpos];
    B_NucValue_old[i] = new double[N_Xpos];
    B_NucValue_tmp[i] = new double[N_Xpos];
    ImmunPopulation[i] = new double[N_Xpos];
    VaccinatedPopulation[i] = new double[N_Xpos];
    SuscepPopRatio[i] = new double[N_Xpos];
    SuscepPopRatio[N_AgeGroup+i] = new double[N_Xpos];    // VaccinatedRatio
    SuscepPopRatio[2*N_AgeGroup+i] = new double[N_Xpos];    // NormalisedVaccinatedRatio
    // AntibodyPopulation
    AntibodyPopulation_All[i] = new double *[N_AllCell];   
     for(j=0; j<N_AllCell; j++)  
     {
      AntibodyPopulation_All[i][j] = new double[jj];
      memset(AntibodyPopulation_All[i][j], 0., jj*SizeOfDouble);
     }
     AntibodyPopulation[i]=AntibodyPopulation_All[i]+DispArray[rank];
   }

  //=============================================================================================================
  //I. Add antibody population - begin
  //============================================================================================================= 
  // First add the antiboday paopulation till previous day infection 
  strcpy(DataFile, "AntigenPopulation.csv");
  for(i=0;i<N_AgeGroup; i++)  // input file has data only for 255 days 
   ReadPopulationData(N_Districts, DataFile, 255, AntibodyPopulation_All[i], InfectRatio[i]);

  // Next add today's antiboday population till previous day infection 
  //Get Int_\Omega_L values using XnodalLdof (assume xnodal & xdof are same)
  if(DistSolOutput) {
   ADISystem3L->IntL(Gamma_Q, IntValue, NucValue, nullptr);}
  else
   {ADISystem3L->IntL(Gamma_Q, IntValue, NucValue, nullptr); }
   
  memcpy(IntValue_old, IntValue, N_Xpos*SizeOfDouble);

  //  if(TDatabase::ParamDB->Par_P0==rank)
  //    for(i=0;i<N_Xpos;i++)    
  //     cout<<" NucValue : "  << NucValue[i]  << " IntValue " << IntValue[i]<< endl;  
  
  // ===================================================================================================
  //testing
  //=====================================================================================================
   TempTotal=0.;

     for(i=0;i<N_Xpos;i++)   
      {
       TempTotal +=IntValue[i];
    //    cout<<" NucValue : "  << NucValue[i]  << " IntValue " << IntValue[i]<< " TempTotal " << TempTotal << endl;
      }

   MPI_Reduce(&TempTotal, &TempTotal_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  
   
   if(TDatabase::ParamDB->Par_P0==rank)         
       cout<<"TempTotal_Actual " << TempTotal_all << " Population ratio " << TempTotal_all/8195.<< endl;

  double *WeibulDist = new double [N_EffAntibdyDays];
  memset(WeibulDist, 0, N_EffAntibdyDays*SizeOfDouble);  
  WeibulDistribution(N_EffAntibdyDays, WeibulDist);

  int today = 0;
  int NewVariant=0;
  // //distribute the initial infected cases across the age group to AntibodyPopulation[]
  // not necessary, day 0 cases will be added after computing B_nuc *******
  // for(i=0;i<N_AgeGroup; i++)
  //  DistributeWaningAntibody(N_Xpos, today, int(end_time), IntValue, WeibulDist, 
  //                           AntibodyPopulation[i], PopRatio[i], nullptr, 1.0, NewVariant);
  // if(TDatabase::ParamDB->Par_P0==rank)
  // for(i=0;i<int(end_time);i++)
  //   cout << i << " WeibulDist[i] " << WeibulDist[i] <<endl; 

  // MPI_Finalize();
  // exit(0);

  //=============================================================================================================
  //I. Add antibody population - end
  //II. Add vaccinated population -begin
  //=============================================================================================================  
  for(i=0;i<N_AgeGroup;i++)
  for(j=0;j<N_Xpos;j++)
   {
    ImmunPopulation[i][j] = AntibodyPopulation[i][j][0];
    VaccinatedPopulation[i][j] = 0.;

    // if(TDatabase::ParamDB->Par_P0==rank)
    // cout << i << ", "<< j<<" ImmunPopulation " <<ImmunPopulation[i][j] <<endl;  
   }
// MPI_Finalize();
//    exit(0);
  // // add vaccined the the respective age group
  // strcpy(DataFile, "InitialVaccinated.dat");
  // ReadData(N_AllCell-2, DataFile, InitialPopulation);

  // //till 23 Mar, only 45+ got vaccinated 
  // //assumed 70% is efefective, (2nd dose considered in data)  
  // for(i=0; i<N_Xpos; i++)
  //  ImmunPopulation[4][i] += InitialPopulation[DispArray[rank] +i]*0.70;      
  
  //=============================================================================================================
  //II. Add vaccinated population - end
  //III. compute Susceptible ratio - begin
  //============================================================================================================= 
  double SL_fact, temp, temp1;
  memset(ImmunPopulation_ALL, 0., N_Xpos*SizeOfDouble);
  for(j=0;j<N_AgeGroup;j++)
   for(i=0; i<N_Xpos; i++)
    ImmunPopulation_ALL[i] += ImmunPopulation[j][i];

  for(j=0;j<N_AgeGroup;j++)
  {
   temp = PopRatio[j];
   for(i=0; i<N_Xpos; i++)
    {
     if(temp>0.) 
     {
      // Antibody
      SuscepPopRatio[j][i] = (temp*StatePopulation[i] - ImmunPopulation[j][i])/(temp*StatePopulation[i]);
      //vaccine
      SuscepPopRatio[N_AgeGroup+j][i] = (temp*StatePopulation[i] - VaccinatedPopulation[j][i])/(temp*StatePopulation[i]); 
     }
     else 
     {
      SuscepPopRatio[j][i] = 0.;
      SuscepPopRatio[N_AgeGroup+j][i] = 0.;
     }
    //  if(TDatabase::ParamDB->Par_P0==rank)
    //   cout << i << ", "<< j<<" SuscepPopRatio " <<ImmunPopulation[j][i] << " " << SuscepPopRatio[j][i] <<endl;   
    } // for i

  //  if(TDatabase::ParamDB->Par_P0==rank)
  //   cout <<endl;

  }// for j

  // MPI_Finalize();
  //  exit(0);
//  if(TDatabase::ParamDB->Par_P0==rank)
//     cout <<endl;

   for(i=0; i<N_Xpos; i++)
   {
  //   temp =0;
  //   for(j=0;j<N_AgeGroup;j++)
  //      temp+= SuscepPopRatio[j][i];

  //  if(TDatabase::ParamDB->Par_P0==rank)
  //   cout<< "temp " << temp <<endl;
   
    // if(temp>0)
    //  SL_fact = ((StatePopulation[i]-ImmunPopulation_ALL[i])/StatePopulation[i])/temp; // S/L
    // else
    //  SL_fact = 1.;

    // for(j=0;j<N_AgeGroup;j++)
    //  if(SuscepPopRatio[j][i]>0)
    //   { SuscepPopRatio[j][i] *= SL_fact; }  
    //  else
    //  { SuscepPopRatio[j][i] = 0.;}

     //normalised VaccinatedRatio
      temp =0;
      for(j=0;j<N_AgeGroup;j++)
       temp+= SuscepPopRatio[N_AgeGroup+j][i];

      for(j=0;j<N_AgeGroup;j++)
       SuscepPopRatio[2*N_AgeGroup+j][i] = double(N_AgeGroup)*SuscepPopRatio[N_AgeGroup+j][i]/temp;

  //  if(TDatabase::ParamDB->Par_P0==rank)
  //   for(j=0;j<N_AgeGroup;j++)
  //     cout << i << ", "<< SL_fact<<" SuscepPopRatio " <<ImmunPopulation[j][i] << " " << SuscepPopRatio[j][i] <<
  //              " Full SuscepPopRatio " <<(StatePopulation[i]-ImmunPopulation_ALL[i])/StatePopulation[i]  <<endl;  

    } // for i

  //  MPI_Finalize();
  //  exit(0);
  //=============================================================================================================
  //III.  compute Susceptible ratio - end
  // Get Nucleation on day 0, distribute to Antibody  - begin
  //============================================================================================================= 

  // ===================================================================================================
  //testing
  //=====================================================================================================
   TempTotal=0.;
   for(i=0;i<N_Xpos;i++)   
    {
     TempTotal +=NucValue[i];
     //  cout<<" NucValue : "  << NucValue[i]  << " IntValue " << IntValue[i]<< " TempTotal " << TempTotal << endl;
    }

   MPI_Reduce(&TempTotal, &TempTotal_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
   if(TDatabase::ParamDB->Par_P0==rank)         
       cout<<"Total_1-Gamma: " << TempTotal_all << endl;
   //  ===================================================================================================

// MPI_Finalize();
//    exit(0);
  Set_sigma_SD(0, 0., 0., 0.000); 

  // NucValue is the actually the \int (1.- \gamma_Q) IdL
  ADISystem3L->Int_B_Nuc(NucValue, B_NucValue_tmp, GrowthAndB_Nuc[0], SuscepPopRatio, N_AgeGroup, MaxAgeOfGroup);


  for(j=0;j<N_AgeGroup;j++)
   {
    memcpy(B_NucValue_old[j], B_NucValue_tmp[j], N_Xpos*SizeOfDouble);
    memset(B_NucValue[j], 0, N_Xpos*SizeOfDouble);  
    
    // add the updated AntibodyPopulation
    for(i=0;i<N_Xpos;i++)
    ImmunPopulation[j][i] -= AntibodyPopulation[j][i][0];
    DistributeWaningAntibody(N_Xpos, today, int(end_time), B_NucValue_tmp[j], WeibulDist, 
                             AntibodyPopulation[j], PopRatio[j], State_CIR, -1.0, NewVariant, StatePopulation);
    for(i=0;i<N_Xpos;i++)
    {
     ImmunPopulation[j][i] += AntibodyPopulation[j][i][0];
    }
   }

#ifdef _MPI

MPI_Request request;
// ==============================================================================
// printing  
// ==============================================================================

// ==============================================================================
// testing 
// ==============================================================================
  memset(TempSum, 0., N_Xpos*SizeOfDouble);
  memset(TempAgeSum, 0., N_AgeGroup*SizeOfDouble);
  
  TempTotal =0;
  for(i=0;i<N_Xpos;i++)
   for(j=0;j<N_AgeGroup;j++)
     {
      TempSum[i] += B_NucValue_tmp[j][i];
      TempAgeSum[j] += B_NucValue_tmp[j][i];
      TempTotal +=B_NucValue_tmp[j][i];
     }

   MPI_Reduce(&TempTotal, &TempTotal_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   

   if(TDatabase::ParamDB->Par_P0==rank)
   for(i=0;i<N_Xpos;i++)
        cout << i<<" B_NucValue_TempSum " << TempSum[i] << " TempTotal_all :" << TempTotal_all <<endl; 

   if(TDatabase::ParamDB->Par_P0==rank)
     cout << "B_NucValueTotal_all :" << TempTotal_all <<endl; 
     
   for(j=0;j<N_AgeGroup;j++)
    {
     MPI_Reduce(TempAgeSum+j, &TempAgeTotal_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

     if(TDatabase::ParamDB->Par_P0==rank)
      cout << i<<" B_AgeRatio " << TempAgeTotal_all/TempTotal_all <<endl;
    }
// ============================================================================================
// Confirmed 
 for(j=0;j<N_AgeGroup;j++)   
  MPI_Gatherv(B_NucValue_tmp[j], N_Xpos, MPI_DOUBLE, CovidPopulation_tmp[j], N_Cells_ALL, DispArray, 
             MPI_DOUBLE, 0, MPI_COMM_WORLD);

 if(rank==0)
  {
   memset(B_NucValue_ALL, 0, N_AllCell*SizeOfDouble);   
   for(j=0;j<N_AgeGroup;j++)    
    {  
     Daxpy(N_AllCell, 1., CovidPopulation_tmp[j], B_NucValue_ALL);  
     WriteInitData(N_AllCell, CovidPopulation_tmp[j], CovidNucleationNameAge[j]); 

    // OutPut("CaseToInfectionRatio " <<  CaseToInfectionRatio <<endl);       
    }

   WriteInitData(N_AllCell, B_NucValue_ALL, CovidNucleationName);                
  }

if(mpi_size>3)
 {
  // Active
  MPI_Gatherv(IntValue, N_Xpos, MPI_DOUBLE, CovidPopulation1, N_Cells_ALL, DispArray, MPI_DOUBLE, 1, MPI_COMM_WORLD);   
  //recovered
  MPI_Gatherv(RecoveredValue, N_Xpos, MPI_DOUBLE, CovidPopulation2, N_Cells_ALL, DispArray, MPI_DOUBLE, 2, MPI_COMM_WORLD);  
  //Deceased
  MPI_Gatherv(RemovedValue, N_Xpos, MPI_DOUBLE, CovidPopulation3, N_Cells_ALL, DispArray, MPI_DOUBLE, 3, MPI_COMM_WORLD);  
  
  if(rank==1)
   {
    WriteInitData(N_AllCell, CovidPopulation1, CovidPopulationName);       
   }
  else if(rank==2)
   {
    WriteInitData(N_AllCell, CovidPopulation2, CovidRecoverdName);    
   }
  else if(rank==3)
   {
    WriteInitData(N_AllCell, CovidPopulation3, CovidDeathName);      
   }
  }
else
  {
   MPI_Gatherv(IntValue, N_Xpos, MPI_DOUBLE, CovidPopulation1, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
   MPI_Gatherv(RecoveredValue, N_Xpos, MPI_DOUBLE, CovidPopulation2, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
   MPI_Gatherv(RemovedValue, N_Xpos, MPI_DOUBLE, CovidPopulation3, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
    
  if(rank==0)
   {
    WriteInitData(N_AllCell, CovidPopulation1, CovidPopulationName);
    WriteInitData(N_AllCell, CovidPopulation2, CovidRecoverdName);    
    WriteInitData(N_AllCell, CovidPopulation3, CovidDeathName);
   }    
 
  // if(rank!=0)
  // {
  //  MPI_Isend(LDistXSum, N_LDistXSum, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);

  // }
  // else
  // {
  //  for(i=1; i<mpi_size; ++i)
  //  {
  //   MPI_Recv(sum_receiv, N_LDistXSum, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  //   for(j=0; j<N_LDistXSum; ++j)
  //    LDistXSum[j] +=sum_receiv[j];
  //  }

  //   WriteInitDistData(N_LnodalPos, LnodalPos, N_LDistXSum, LDistXSum);  
  // }

  // if(rank==0)
  //  {
  //   WriteInitRecovDistData(N_LnodalPos, LnodalPos, N_LDistXSum, LRecoDistXSum);  
  //  }
  }
  //  OutPut("CaseToInfectionRatio " <<  CaseToInfectionRatio <<endl);  
  //  MPI_Finalize();
  //  exit(0);
#else
//   // for(i=0; i<N_Xpos; i++)
//   //  OutPut("Covid Population in " <<  enum_to_string(i) << " : " << IntValue[i] <<endl);
//  WriteInitData(N_Xpos, IntValue, CovidPopulationName);
//  WriteInitData(N_Xpos, B_NucValue_ALL, CovidNucleationName);
//   // //  for(j=0;j<N_AgeGroup;j++)              
//  WriteInitData(N_Xpos, B_NucValue[0], CovidNucleationName0);
//  WriteInitData(N_Xpos, B_NucValue[1], CovidNucleationName1);
//  WriteInitData(N_Xpos, B_NucValue[2], CovidNucleationName2);  
//  WriteInitData(N_Xpos, RecoveredValue, CovidRecoverdName);
//  WriteInitData(N_Xpos, RemovedValue, CovidDeathName);
#endif

  //============================================================================================================= 
  //III. compute Susceptible ratio before time loop - begin
  //============================================================================================================= 
  memset(ImmunPopulation_ALL, 0., N_Xpos*SizeOfDouble);
  AntibodyPop = 0;  
  VaccinePop = 0;
  for(j=0;j<N_AgeGroup;j++)
   for(i=0; i<N_Xpos; i++)
    {
     ImmunPopulation_ALL[i] += ImmunPopulation[j][i];
     AntibodyPop += ImmunPopulation[j][i];
     VaccinePop += VaccinatedPopulation[j][i];
    }

   MPI_Reduce(&AntibodyPop, &AntibodyPopTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
   MPI_Reduce(&VaccinePop, &VaccinePopTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  

   if(TDatabase::ParamDB->Par_P0==rank)
    {
     WriteInitAntiBody(AntibodyPopTot, "Antibody");
     WriteInitAntiBody(VaccinePopTot, "Vaccinated");
     OutPut("Day : "<< 0 <<"AntibodyPopTot "  << AntibodyPopTot <<" VaccinePopTot "  << VaccinePopTot <<endl);
    }

  for(j=0;j<N_AgeGroup;j++)
   {
    temp = PopRatio[j];
   for(i=0; i<N_Xpos; i++)
    {
     if(temp>0.) 
     { 
       //Antibody
      SuscepPopRatio[j][i] = (temp*StatePopulation[i] - ImmunPopulation[j][i])/(temp*StatePopulation[i]);
      //vaccine
      SuscepPopRatio[N_AgeGroup+j][i] = (temp*StatePopulation[i] - VaccinatedPopulation[j][i])/(temp*StatePopulation[i]); 
     }
     else 
     {
      SuscepPopRatio[j][i] = 0.;
      SuscepPopRatio[N_AgeGroup+j][i] =0.; 
     }
    } // for i
   }// for j

  for(i=0; i<N_Xpos; i++)
   {
    // temp =0;
    // for(j=0;j<N_AgeGroup;j++)
    //    temp+= SuscepPopRatio[j][i];

    // if(temp>0)
    //  SL_fact = ((StatePopulation[i]-ImmunPopulation_ALL[i])/StatePopulation[i])/temp; // S/L
    // else
    //  SL_fact = 1.;

    // for(j=0;j<N_AgeGroup;j++)
    //  if(SuscepPopRatio[j][i]>0)
    //   { SuscepPopRatio[j][i] *= SL_fact; }  
    //  else
    //  { SuscepPopRatio[j][i] = 0.;}

     //normalised VaccinatedRatio
      temp =0;
      for(j=0;j<N_AgeGroup;j++)
       temp+= SuscepPopRatio[N_AgeGroup+j][i];

      for(j=0;j<N_AgeGroup;j++)
       SuscepPopRatio[2*N_AgeGroup+j][i] = double(N_AgeGroup)*SuscepPopRatio[N_AgeGroup+j][i]/temp;     
    } // for i
  //=============================================================================================================
  //III.  compute Susceptible ratio - end
  //=============================================================================================================
  

  // //  MPI_Finalize();
  // //  exit(0);

  
#ifdef _MPI
 t1 = MPI_Wtime(); 
#endif

  double NoNucTime = 140.5; // first sun locdown day since Mar 23
  // double NoNucIncr=3.;
  //======================================================================
  // time disc loop
  //======================================================================    
   UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
   UpdateRhs = TRUE; //check BilinearCoeffs in example file
   ConvectionFirstTime=TRUE;

   // time loop starts
   while(TDatabase::TimeDB->CURRENTTIME < end_time)
   {
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(l=0; l<N_SubSteps; l++) // sub steps of fractional step theta
    {
      SetTimeDiscParameters(1);

#ifdef _MPI
 if(TDatabase::ParamDB->Par_P0==rank)
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
 
      //  MPI_Finalize();
      //  exit(0);
      //set parameters for nucleation
      // after 93 days sd value increases by 0.02 for every 15 days,
      Set_sigma_SD(TDatabase::TimeDB->CURRENTTIME, NoNucTime, 15., 0.000);     
      // if(TDatabase::TimeDB->CURRENTTIME-0.5>NoNucTime)
      //  { 
      //   NoNucTime +=NoNucIncr;
      //   if(NoNucIncr==3.) // no nuc on Sun & Wed 
      //    {NoNucIncr=4.;}
      //   else 
      //   {NoNucIncr=3.;}
      //  } 
// #ifdef _MPI 
//  if(TDatabase::ParamDB->Par_P0==rank)
//  #endif        
//     OutPut(endl << "CURRENT TIME: " << TDatabase::TimeDB->CURRENTTIME << endl); 
    // if(rank==0)
      // OutPut("NoNucTime: " <<NoNucTime<< endl);    
    //  cout << "IntVal L start "<<endl;
    //  //Get Int_\Omega_L values using XnodalLdof (assume xnodal & xdof are same)
    //  ADISystem3L->IntL(Gamma_Q, IntValue, NucValue);
    //  cout << "IntVal L end "<<endl;
     // XnodalLnodal->XnodalLdof->LdofXdof (note XnodalLdof is needed for IntL)

      // time-dep reaction
      GetReactionFactor(N_LnodalPos, LnodalPos, ReactionFactor);

      ADISystem3L->XnodalLnodalToLdofXdof(MaxN_PtsForNodal, Sol_LdofXdof);
       
      #ifdef __LOCKDOWNMODEL__
       for(i=0; i<N_LNodals_All; ++i)
        {
         Sol_Xdof = Sol_LdofXdof + i*N_DOF;
         Recovered_Xdof =  Recovery_LdofXdof + i*N_DOF;
         GrowthFact_X = (1. - tau*TDatabase::TimeDB->THETA2*ReactionFactor[2*i])/(1. + tau*TDatabase::TimeDB->THETA1*ReactionFactor[2*i]);   
         RecoveryRatio = ReactionFactor[2*i+1]/ReactionFactor[2*i]; //recovery ratio in the total reaction

         for(jj=0; jj<N_DOF; ++jj)
          {            
            // I_R = RecoveryRatio * Totalremoved = RecoveryRatio*(I^n - I^{n+1})
           Recovered_Xdof[jj] = RecoveryRatio*(1.-GrowthFact_X)*Sol_Xdof[jj]; // 
           Sol_Xdof[jj] = GrowthFact_X*Sol_Xdof[jj];        
          }
        }
      #else
        cout <<"NON __LOCKDOWNMODEL__ Not yet implemented " <<endl;
        exit(0);
      #endif
 
    //===============================================================================
    // compute the recovery values
    //===============================================================================
     // transfer/copy reovery values from Recovery_LdofXdof to Recovery_XnodalLdof
     ADISystem3L->LdofXdofToXnodalLdof(Recovery_LdofXdof); 

     //Get Int_\Omega_L values using XnodalLdof (assume xnodal & xdof are same)
     //NucValue value (3rd argument) will not be computed when first argument is NULL
     //Recoverd solution is in the system, so Int_L will compute the total recoverd
     ADISystem3L->IntL(nullptr, RecoveredValue_temp, NucValue, nullptr); 
     Daxpy(N_Xpos, 1., RecoveredValue_temp, RecoveredValue); 
    //  Daxpy(N_LDistXSum, 1., LRecoDistXSum_temp, LRecoDistXSum);          
    //===============================================================================
  
     // transfer sol from Sol_LdofXdof to Sol_XnodalLdof
     ADISystem3L->LdofXdofToXnodalLdof(Sol_LdofXdof); 
     //===============================================================================
     //system X-direction solution -- end
     //system L-direction solution -- start  
     //===============================================================================
  
     //Get I and (1.- \gamma_Q)I over Int_\Omega_L using XnodalLdof (assume xnodal & xdof are same)
     ADISystem3L->IntL(Gamma_Q, IntValue, NucValue, NULL);

     //compute death cases
     Daxpy(N_Xpos, -1., RecoveredValue_temp, IntValue_old); //   y := alpha*x + y 
     DsumS(N_Xpos, 1., IntValue_old, IntValue, RemovedValue);// Infect Death 
     
     //===============================================================================
     //Solve L0-direction solution -- start  
     //===============================================================================
#ifdef __LD__     
  // // ===================================================================================================
  // //testing
  // //=====================================================================================================
  //  TempTotal=0.;
  //  for(i=0;i<N_Xpos;i++)   
  //   {
  //    TempTotal +=NucValue[i];

  //    if(rank==4)
  //    printf( "Rank %d NucValue :%f IntValue : %f TempTotal: %f \n", rank, NucValue[i], IntValue[i], TempTotal );    
  //   }

  // // printf("Rank %d TempTotal: %f \n", rank,  TempTotal );

  //  MPI_Reduce(&TempTotal, &TempTotal_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
  //  if(TDatabase::ParamDB->Par_P0==rank)         
  //      cout<<"Total_1-Gamma: " << TempTotal_all << endl;
  // //===================================================================================================
  // //testing       
  // //===================================================================================================
      
      // NucValue is the actually the \int (1.- \gamma_Q)I   
      // compute B_Nuc only for output purpose              
      ADISystem3L->Int_B_Nuc(NucValue, B_NucValue_tmp, GrowthAndB_Nuc[0], SuscepPopRatio, N_AgeGroup, MaxAgeOfGroup);
  
      for(j=0;j<N_AgeGroup;j++)
       {
        DsumP(N_Xpos, 0.5*tau, B_NucValue_old[j], B_NucValue_tmp[j], B_NucValue[j]);
        memcpy(B_NucValue_old[j], B_NucValue_tmp[j], N_Xpos*SizeOfDouble);
      }

// // ==============================================================================
// // testing 
// // ==============================================================================
//   memset(TempSum, 0., N_Xpos*SizeOfDouble);
//   memset(TempAgeSum, 0., N_AgeGroup*SizeOfDouble);
  
//   TempTotal =0;
//   for(i=0;i<N_Xpos;i++)
//    for(j=0;j<N_AgeGroup;j++)
//      {
//       // TempSum[i] += B_NucValue_tmp[j][i];
//       // TempAgeSum[j] += B_NucValue_tmp[j][i];
//       // TempTotal +=B_NucValue_tmp[j][i];
//       TempSum[i] += B_NucValue[j][i];
//       TempAgeSum[j] += B_NucValue[j][i];
//       TempTotal +=B_NucValue[j][i];      
//      }

//    MPI_Reduce(&TempTotal, &TempTotal_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   

//   //  if(TDatabase::ParamDB->Par_P0==rank)
//   //  for(i=0;i<N_Xpos;i++)
//         // cout << i<<" B_NucValue_TempSum " << TempSum[i] << " TempTotal_all :" << TempTotal_all <<endl; 

//   //  if(TDatabase::ParamDB->Par_P0==rank)
//   //    cout << "B_NucValueTotal_all :" << TempTotal_all <<endl; 
     
//    for(j=0;j<N_AgeGroup;j++)
//     {
//      MPI_Reduce(TempAgeSum+j, &TempAgeTotal_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//      if(TDatabase::ParamDB->Par_P0==rank)
//       cout << i<< " Age B_Nuc: " << TempAgeTotal_all << " B_AgeRatio " << TempAgeTotal_all/TempTotal_all <<endl;
//     }

//    if(TDatabase::ParamDB->Par_P0==rank)
//     cout <<endl;
//  if(m==1)
//    { MPI_Finalize();
//    exit(0); }    
// // ===================================================================================================
// //testing  
// // ============================================================================================


     // ADISystem3L->CopySolToInternal(0); // do nothing for ADISystem3L[0]
     //solve inernal system for every x pos
     // B_Nuc will actually be included in the model
     ADISystem3L->Solve(0, BilinearCoeffs_L0, NucValue, SuscepPopRatio, N_AgeGroup, MaxAgeOfGroup);

    //  ADISystem3L->CopySolFromInternal(0); // do nothing for ADISystem3L[0]     

// //=========================================================================
// //testing
// //=========================================================================
//   //only second argument will be calculated
//   ADISystem3L->IntL(nullptr, IntValue, NucValue, nullptr);

//    TempTotal=0.;
//    for(i=0;i<N_Xpos;i++)   
//     {
//      TempTotal +=IntValue[i];
//      //  cout<<" NucValue : "  << NucValue[i]  << " IntValue " << IntValue[i]<< " TempTotal " << TempTotal << endl;
//     }

//    MPI_Reduce(&TempTotal, &TempTotal_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
//    if(TDatabase::ParamDB->Par_P0==rank)         
//        cout<<"I val_all: " << TempTotal_all << endl;
// // if(m==2)
// //    { MPI_Finalize();
// //    exit(0); }        
// //=========================================================================
// //testing
// //=========================================================================

#endif    
 
     //===============================================================================
     //Solve L0-direction solution -- completed  
     //Solve L1-direction solution -- start  
     //===============================================================================
#ifdef __LV__  
     // XnodalL2L1L0nodal -->  XnodalL2L0L1nodal       
     ADISystem3L->CopySolToInternal(1);
   
     // solve inernal system for every x pos
     ADISystem3L->Solve(1, BilinearCoeffs_L1, IntValue, SuscepPopRatio, N_AgeGroup, MaxAgeOfGroup);

     // XnodalL2L0L1nodal --> XnodalL2L1L0nodal
     ADISystem3L->CopySolFromInternal(1);  
  // OutPut(" __LV__ Complete "  << endl);       
#endif        
     //===============================================================================
     //Solve L1-direction solution -- completed  
     //Solve L2-direction solution -- start  
     //===============================================================================
#ifdef __LA__      
     // XnodalL2L1L0nodal  --> XnodalL0L1L2nodal
     ADISystem3L->CopySolToInternal(2);

     // solve inernal system for every x pos
     ADISystem3L->Solve(2, BilinearCoeffs_L2, IntValue, SuscepPopRatio, N_AgeGroup, MaxAgeOfGroup);
    
     // XnodalL0L1L2nodal --> XnodalL2L1L0nodal
     ADISystem3L->CopySolFromInternal(2);
  // OutPut(" __LA__ Complete "  << endl);          
#endif     
     //===============================================================================
     //Solve L2-direction solution -- completed  
     //===============================================================================
    } // for(l=0;l<N_SubSteps;l++) 
    
     //Get Int_\Omega_L values using XnodalLdof (assume xnodal & xdof are same)   
     if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0 && DistSolOutput)   
     {
      ADISystem3L->IntL(Gamma_Q, IntValue, NucValue, nullptr);
     }
     else
     {
      ADISystem3L->IntL(Gamma_Q, IntValue, NucValue, nullptr);
     }

     memcpy(IntValue_old, IntValue, N_Xpos*SizeOfDouble);

    //======================================================================
    // produce outout
    //======================================================================
    if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
    {
     #ifdef _MPI

// // ==========================================================================
// // testing-begin
// // ==========================================================================
//      TempTotal =0;
//      memset(TempAgeSum, 0., N_AgeGroup*SizeOfDouble);
//      for(i=0;i<N_Xpos;i++)
//       for(j=0;j<N_AgeGroup;j++)
//        {
//         // TempSum[i] += B_NucValue_tmp[j][i];
//         TempAgeSum[j] += B_NucValue[j][i];
//         // TempTotal +=B_NucValue_tmp[j][i];
//        }
         
//    for(j=0;j<N_AgeGroup;j++)
//     {
//      MPI_Reduce(TempAgeSum+j, &TempAgeTotal_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//      if(TDatabase::ParamDB->Par_P0==rank)
//       cout << i<< " Age B_Nuc II: " << TempAgeTotal_all << " B_AgeRatio " << TempAgeTotal_all/TempTotal_all <<endl;
//     }
//    if(TDatabase::ParamDB->Par_P0==rank)
//     cout <<endl;
// // ==========================================================================
// // testing-end
// // ==========================================================================

      //confirmed
      for(j=0;j<N_AgeGroup;j++)   
        MPI_Gatherv(B_NucValue[j], N_Xpos, MPI_DOUBLE, CovidPopulation_tmp[j], N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if(rank==0)
       {
        memset(B_NucValue_ALL, 0, N_AllCell*SizeOfDouble);   
         for(j=0;j<N_AgeGroup;j++)   
          { 
           Daxpy(N_AllCell, 1., CovidPopulation_tmp[j], B_NucValue_ALL);  
           WriteData(N_AllCell, CovidPopulation_tmp[j], CovidNucleationNameAge[j]); 

          //  TempTotal=0.;
          //   for(i=0;i<N_AllCell;i++)   
          //     TempTotal +=CovidPopulation_tmp[j][i];

          //    cout<<"B_NucAgeGroup_All : " << TempTotal << " B_AgeRatio " << TempTotal/TempTotal_all << endl;
          } 
          // cout <<endl;
         WriteData(N_AllCell, B_NucValue_ALL, CovidNucleationName);                
       }  
 
      if(mpi_size>3)
       {
        MPI_Gatherv(IntValue, N_Xpos, MPI_DOUBLE, CovidPopulation1, N_Cells_ALL, DispArray, MPI_DOUBLE, 1, MPI_COMM_WORLD);   
        MPI_Gatherv(RecoveredValue, N_Xpos, MPI_DOUBLE, CovidPopulation2, N_Cells_ALL, DispArray, MPI_DOUBLE, 2, MPI_COMM_WORLD);  
        MPI_Gatherv(RemovedValue, N_Xpos, MPI_DOUBLE, CovidPopulation3, N_Cells_ALL, DispArray, MPI_DOUBLE, 3, MPI_COMM_WORLD);  
  
        if(rank==1) { 
          //  TempTotal=0.;
          //   for(i=0;i<N_AllCell;i++)   
          //     TempTotal +=CovidPopulation1[i];

          //    cout<<"CovidPopulation1 : " << TempTotal << endl;   

          WriteData(N_AllCell, CovidPopulation1, CovidPopulationName);  }
        else if(rank==2) { 
          WriteData(N_AllCell, CovidPopulation2, CovidRecoverdName); }
        else if(rank==3) {
          WriteData(N_AllCell, CovidPopulation3, CovidDeathName); }
       }
      else {
        MPI_Gatherv(IntValue, N_Xpos, MPI_DOUBLE, CovidPopulation1, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
        MPI_Gatherv(RecoveredValue, N_Xpos, MPI_DOUBLE, CovidPopulation2, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
        MPI_Gatherv(RemovedValue, N_Xpos, MPI_DOUBLE, CovidPopulation3, N_Cells_ALL, DispArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);  

        if(rank==0) {
          WriteData(N_AllCell, CovidPopulation1, CovidPopulationName);
          WriteData(N_AllCell, CovidPopulation2, CovidRecoverdName);    
          WriteData(N_AllCell, CovidPopulation3, CovidDeathName);
        }     
      }
     #else
      // WriteData(N_Xpos, IntValue, CovidPopulationName);   
      // WriteData(N_Xpos, B_NucValue[0], CovidNucleationName0);  
      // WriteData(N_Xpos, B_NucValue[1], CovidNucleationName1);  
      // WriteData(N_Xpos, B_NucValue[2], CovidNucleationName2);              
      // WriteData(N_Xpos, RecoveredValue, CovidRecoverdName);           
      // WriteData(N_Xpos, RemovedValue, CovidDeathName);     
     #endif

     today++;
     if(NewVariant>0)
      NewVariant++;

     for(j=0;j<N_AgeGroup;j++)      
      DistributeWaningAntibody(N_Xpos, today, int(end_time), B_NucValue[j], WeibulDist, 
                               AntibodyPopulation[j], PopRatio[j], State_CIR, -1.0, NewVariant, StatePopulation);
 
  //============================================================================================================= 
  // vaccine schedule  - begin  
  //=============================================================================================================      
     // 1 Aug - 31 DEC = 153 days , 245
     // if(TDatabase::TimeDB->CURRENTTIME>245) {
       // IND: 0.0027=33L per day and 70% 
    //============================================================================================================= 
    // KA Vaccination policy as on 26th Jun 2,13,42,378
    //============================================================================================================= 
       //Phase I, age>45, till 4 Apr, KA: 0.00363006341  = (4608184/(0.2325*70000000))/78 
       //Phase II, from 5th Apr: 0.00423253479 = (( 21342378- 4608184)/((0.2325+0.448)*70000000))/83
       //KA: 0.0019 = (21342378/70000000)/160 till  24th Jun ==> 1,33,389 doses per day
       //As on 26th Jun: 21342378= 0.2325*70000000*0.00363006341*78 + (0.2325+0.448)*70000000*0.00423253479*83

       //KA : 0.00442396531 = 9117535/(0.448*70774100*65) // 18 - 44 yrs, from 1 May to 4 Jul (65 days)
       //KA : 0.0092558902 = 8158503/(0.136860*70774100*91) // 45 - 59 yrs, from 4 Apr to 4 Jul (91 days)       
       //KA : 0.00764903233 = 6516629/(0.095537*70774100*126) // Above 60 yrs, from 1 Mar to 4 Jul (31 days)
    //=============================================================================================================        
    if(today<365)
      { AcclerateVacRate = 1.; }
    else
      { AcclerateVacRate = TDatabase::ParamDB->KDP_c_p; }

      // 18 - 44 yrs, from 1 May (304 days since 1 Jul 2020)     
      if(TDatabase::TimeDB->CURRENTTIME>304)  
        for(i=0; i<N_Xpos; i++) 
         {
          VaccinatedPopulation[2][i] +=  PopRatio[2]*StatePopulation[i]*0.00442396531*0.35*AcclerateVacRate; // 18<age<45
          if(VaccinatedPopulation[2][i] > PopRatio[2]*StatePopulation[i])
             VaccinatedPopulation[2][i] = PopRatio[2]*StatePopulation[i]; // 100% vaccinated
         }

     // 45 - 59 yrs, from 4 Apr  (277 days since 1 Jul 2020)   
     if(TDatabase::TimeDB->CURRENTTIME>277)  
        for(i=0; i<N_Xpos; i++)          
         {
          VaccinatedPopulation[3][i] += PopRatio[3]*StatePopulation[i]*0.0092558902*0.35*AcclerateVacRate; // age>45       
           if( VaccinatedPopulation[3][i]> PopRatio[3]*StatePopulation[i]) 
              VaccinatedPopulation[3][i] = PopRatio[3]*StatePopulation[i]; // 100% vaccinated
         }
     
     //From 1 Mar 2021, // above 60 (243 days since 1 Jul 2020)  
     if(TDatabase::TimeDB->CURRENTTIME>243)     
       for(i=0; i<N_Xpos; i++)          
        {     
         VaccinatedPopulation[4][i] += PopRatio[4]*StatePopulation[i]*0.00764903233*0.35*AcclerateVacRate; // age>45       
          if( VaccinatedPopulation[4][i]> PopRatio[4]*StatePopulation[i]) 
             VaccinatedPopulation[4][i] = PopRatio[4]*StatePopulation[i]; // 100% vaccinated
        }
  //============================================================================================================= 
  // vaccine schedule  - end
  //introduce new variant
  //=============================================================================================================   
   if(NewVariant_day==today || today==NewVariant_day_Wave3)
    {
     InitNewVirusVariant(today, N_AgeGroup, N_Xpos, N_EffAntibdyDays, AntibodyPopulation);
     NewVariant =1;     
    }
  //============================================================================================================= 
    // Compute Today's AB - begin
  //============================================================================================================= 
    for(j=0;j<N_AgeGroup;j++)   
     for(i=0; i<N_Xpos; i++) 
      {
       temp = PopRatio[j]*StatePopulation[i];
       temp1 = VaccinatedPopulation[j][i] + AntibodyPopulation[j][i][today];
       
       if(temp1<temp) {
         ImmunPopulation[j][i] = temp1; }
       else {
         ImmunPopulation[j][i] = temp; }
      }
  //============================================================================================================= 
  // Compute Today's AB - end  
  // compute Susceptible ratio - begin
  //============================================================================================================= 
     memset(ImmunPopulation_ALL, 0., N_Xpos*SizeOfDouble);
     AntibodyPop = 0; 
     VaccinePop = 0; 
     for(j=0;j<N_AgeGroup;j++)
     for(i=0; i<N_Xpos; i++) 
     {
      ImmunPopulation_ALL[i] += ImmunPopulation[j][i];
      AntibodyPop += ImmunPopulation[j][i]; 
      VaccinePop += VaccinatedPopulation[j][i];
      }

     MPI_Reduce(&AntibodyPop, &AntibodyPopTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   
     MPI_Reduce(&VaccinePop, &VaccinePopTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

     if(rank==0)
      {
       WriteAntiBody(AntibodyPopTot, "Antibody");
       WriteAntiBody(VaccinePopTot, "Vaccinated");
       //OutPut("Day : "<< today <<" AntibodyPopTot "  << AntibodyPopTot <<" VaccinePopTot "  << VaccinePopTot <<endl);
      }

     for(j=0;j<N_AgeGroup;j++) {
        temp = PopRatio[j];
      for(i=0; i<N_Xpos; i++) {
        if(temp>0.) {
          SuscepPopRatio[j][i] = (temp*StatePopulation[i] - ImmunPopulation[j][i])/(temp*StatePopulation[i]);
          SuscepPopRatio[N_AgeGroup+j][i] = (temp*StatePopulation[i] - VaccinatedPopulation[j][i])/(temp*StatePopulation[i]); 
         }
        else  
        {
          SuscepPopRatio[j][i] = 0.;
          SuscepPopRatio[N_AgeGroup+j][i] = 0.; 
        }
       } // for i 
      }// for j

    for(i=0; i<N_Xpos; i++) 
     {
      // temp =0;
      // for(j=0;j<N_AgeGroup;j++)
      //  temp+= SuscepPopRatio[j][i];

      // if(temp>0)
      //   SL_fact = ((StatePopulation[i]-ImmunPopulation_ALL[i])/StatePopulation[i])/temp; // S/L
      // else
      //   SL_fact = 1.;

      // for(j=0;j<N_AgeGroup;j++)
      //   if(SuscepPopRatio[j][i]>0) { 
      //    SuscepPopRatio[j][i] *= SL_fact; }  
      //   else { 
      //    SuscepPopRatio[j][i] = 0.; }


     //normalised VaccinatedRatio
      temp =0;
      for(j=0;j<N_AgeGroup;j++)
       temp+= SuscepPopRatio[N_AgeGroup+j][i];

      for(j=0;j<N_AgeGroup;j++)
       SuscepPopRatio[2*N_AgeGroup+j][i] = double(N_AgeGroup)*SuscepPopRatio[N_AgeGroup+j][i]/temp;         
     } // for i
  //=============================================================================================================
  //III.  compute Susceptible ratio - end
  // reset Bnuc and other daily data
  //=============================================================================================================
 
    for(i=0; i<N_AgeGroup; i++) {
      memset(B_NucValue[i], 0, N_Xpos*SizeOfDouble);  
     }
           
     memset(RemovedValue, 0, N_Xpos*SizeOfDouble);  
     memset(RecoveredValue, 0, N_Xpos*SizeOfDouble);     
     memset(LRecoDistXSum, 0, N_LDistXSum*SizeOfDouble);   
    }   // if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
   } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

#ifdef _MPI
 t2 = MPI_Wtime();
 if(TDatabase::ParamDB->Par_P0==rank)
  printf( "Elapsed time: %f\n", t2 - t1 ); 
#endif  

CloseFiles();
  
#ifdef _MPI
MPI_Finalize();
#endif

return 0;
} // end main
