/** ************************************************************************ 
* @brief     source file for TSystemTCD3D_3DSuraface
* @author    Sashikumaar Ganesan
* @date      07.10.2019
* @History 
 ************************************************************************  */
 #include <Database.h>
 #include <SystemRTE.h>
 #include <SquareStructure3D.h>
 #include <DiscreteForm3D.h>
 #include <Assemble3D.h>
 #include <AuxParam3D.h>
 #include <MultiGridScaIte.h>
 #include <LocalProjection.h>
 #include <DirectSolver.h>
 #include <Solver.h>
 #include <AssembleMat3D.h>
 #include <FEVectFunct3D.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
 #include <FEDatabase3D.h>

#include <stdlib.h>
#include <string.h>

TSystemRTE::TSystemRTE(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver,
                       TFEFunction3D **spatialFEFunctions, DoubleFunctND *getKernel, DoubleFunctND *getRhs)
                       :TSystemTCD3D(N_levels, fespaces, sol, rhs, disctype, solver)
{
 char  *GEO, *PRM;
 int N_Cells, N_dof_level;
#ifdef _SMPI
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  MPI_Comm_size(TDatabase::ParamDB->Comm, &mpi_size);  
#endif 
 //=========================================================================
 // Internal space
 // construct FESpace and matrices for internal space matrices
 //=========================================================================
 
  PRM = TDatabase::ParamDB->BNDFILE_INTL;
  GEO = TDatabase::ParamDB->GEOFILE_INTL;
  char *SMESH = TDatabase::ParamDB->SMESHFILE;

  //generate mesh
   Domain_Intl = new TDomain();   
  if (!strcmp(GEO, "TetGenMesh"))
  { 
   Domain_Intl->TetrameshGen(PRM, SMESH);  // Use Tetgen 
  }
  else
  {
   Domain_Intl->Init(PRM, GEO); // ParMooN  build-in Geo mesh
  }

  LEVELS_INTL = TDatabase::ParamDB->LEVELS_INTL;
  if(TDatabase::ParamDB->SOLVER_TYPE_INTL==DIRECT)
  {
    TDatabase::ParamDB->UNIFORM_STEPS_INTL += (LEVELS_INTL-1);
    LEVELS_INTL = 1;
  }

  for(int i=0;i<TDatabase::ParamDB->UNIFORM_STEPS_INTL;i++)
    Domain_Intl->RegRefineAll();
  
  Coll_Intl=Domain_Intl->GetCollection(It_Finest, 0);
  N_Cells = Coll_Intl->GetN_Cells();
  Cell_IDX = new int[N_Cells];
  Joint_IDX = new int[6*N_Cells];

  //generate a surface 2d mesh using Domain_Intl
  SurfDomain = Domain_Intl->GetSurfMesh(Cell_IDX, Joint_IDX);
  Surf_Coll = SurfDomain->GetCollection(It_Finest, 0);
  N_SurfCells = Surf_Coll->GetN_Cells();

 #ifdef _SMPI
 int i;
  for(i=0;i<N_SurfCells;++i)
   Surf_Coll->GetCell(i)->SetGlobalCellNo(i);
     
  N_AllSurfCell = N_SurfCells;


 if(N_SurfCells<mpi_size)
  {
   if(TDatabase::ParamDB->Par_P0==1)
    printf("N_SurfCells less than N_MPI_Processes: N_Cells = %d\n", N_Cells);
   MPI_Finalize();
   exit(0); 
  }
  
  int N_OwnSurfCells = N_SurfCells/mpi_size; 
  int N_RemainingCells = N_SurfCells%mpi_size;   

  if(TDatabase::ParamDB->Par_P0==1)
  printf("N_OwnSurfCells = %d, N_RemainingCells = %d\n", N_OwnSurfCells, N_RemainingCells);

  int AdditionCell=0;
  if(rank<N_RemainingCells)
   AdditionCell = 1 ;

  N_SurfCells = N_OwnSurfCells + AdditionCell;

  TBaseCell **OwnCells;
  OwnCells = new TBaseCell*[N_SurfCells];
  int *GLOB_cellIndex = new int[N_SurfCells];

  N_Cells_ALL = new int[mpi_size];

  MPI_Allgather(&N_SurfCells, 1, MPI_INT, N_Cells_ALL, 1, MPI_INT, TDatabase::ParamDB->Comm);

  int disp = 0;
  for(i=0;i<rank;++i)
   disp +=N_Cells_ALL[i];

  for(i=0;i<N_SurfCells;++i)
  {
   OwnCells[i] = Surf_Coll->GetCell(disp+i);
   GLOB_cellIndex[i] = disp+i;
  }

  SurfDomain->ReplaceTreeInfo(N_SurfCells, OwnCells, GLOB_cellIndex, N_SurfCells);
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);

  Surf_Coll = SurfDomain->GetCollection(It_Finest, 0);
  N_SurfCells = Surf_Coll->GetN_Cells();
 #endif

  // cout<< "N_SurfCells " << N_SurfCells <<endl;

  B1 = new TSquareMatrix3D*[N_Levels];
  B2 = new TSquareMatrix3D*[N_Levels];
  B3 = new TSquareMatrix3D*[N_Levels];

  if (Disctype==SUPG)
   {
     S1_SUPG = new TSquareMatrix3D*[N_Levels];
     S2_SUPG = new TSquareMatrix3D*[N_Levels];
     S3_SUPG = new TSquareMatrix3D*[N_Levels];
     S1S2_SUPG = new TSquareMatrix3D*[N_Levels];
     S1S3_SUPG = new TSquareMatrix3D*[N_Levels];
     S2S3_SUPG = new TSquareMatrix3D*[N_Levels];
     M1_SUPG = new TSquareMatrix3D*[N_Levels];     
     M2_SUPG = new TSquareMatrix3D*[N_Levels];
     M3_SUPG = new TSquareMatrix3D*[N_Levels];
   }
 

  for(int i=Start_Level;i<N_Levels;i++)
  {
   B1[i] = new TSquareMatrix3D(sqstructure[i]);
   B2[i] = new TSquareMatrix3D(sqstructure[i]); 
   B3[i] = new TSquareMatrix3D(sqstructure[i]); 

   if (Disctype==SUPG)
   {
     S1_SUPG[i] = new TSquareMatrix3D(sqstructure[i]);
     S2_SUPG[i] = new TSquareMatrix3D(sqstructure[i]);
     S3_SUPG[i] = new TSquareMatrix3D(sqstructure[i]);
     S1S2_SUPG[i] = new TSquareMatrix3D(sqstructure[i]);
     S1S3_SUPG[i] = new TSquareMatrix3D(sqstructure[i]);
     S2S3_SUPG[i] = new TSquareMatrix3D(sqstructure[i]);
     M1_SUPG[i] = new TSquareMatrix3D(sqstructure[i]);
     M2_SUPG[i] = new TSquareMatrix3D(sqstructure[i]);
     M3_SUPG[i] = new TSquareMatrix3D(sqstructure[i]);          
   }

  }

  Rhs1Array = new double *[N_Levels];
  Rhs2Array = new double *[N_Levels];
  Rhs3Array = new double *[N_Levels];

   if (Disctype==SUPG)
   for(int i=0;i<9;i++)
    Rhs_SUPG[i] = new double *[N_Levels];

  for(int i=Start_Level;i<N_Levels;i++)
  {
   N_dof_level=FeSpaces[i]->GetN_DegreesOfFreedom();
   Rhs1Array[i] = new double[N_dof_level];   
   Rhs2Array[i] = new double[N_dof_level];    
   Rhs3Array[i] = new double[N_dof_level];   
  }

   if (Disctype==SUPG)
   {
    for(int j=0;j<9;j++)
     for(int i=Start_Level;i<N_Levels;i++)
     Rhs_SUPG[j][i]  = new double[N_dof_level];   
   }

  rhs1old =  new double[N_DOF];
  rhs2old =  new double[N_DOF];
  rhs3old =  new double[N_DOF];

  RhsAssemble = new TAssembleMat3D*[N_Levels];

  /** need it for solver */
  sqmatrices = (TSquareMatrix **)SQMATRICES_RTE;

  SpatialFEFunctions = spatialFEFunctions;

  GetKernel = getKernel;
  GetRhs = getRhs;
  // cout << "Domain_Intl " << N_SurfCells <<endl;

  // exit(0);
  
} // constructor


TSystemRTE::~TSystemRTE()
{
  int i;
  
  for(i=Start_Level;i<N_Levels;i++)
   {
    delete sqstructure[i];
    delete sqmatrixA[i];   
   }
   
    delete [] sqstructure;
    delete [] sqmatrixA;
  
  if (SOLVER==GMG && TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
   {
    delete [] Itmethod_sol;
    delete [] Itmethod_rhs;
   }
  
  delete [] B;
  delete [] defect;
  
//   delete sqmatrixM;
//   if(Disctype==SDFEM || Disctype==SUPG)
//    {
//      delete sqmatrixS;
//      delete sqmatrixK;    
//    }
}


 void TSystemRTE::Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue,
                            BoundCondFunct3D *BoundCond_Intl, BoundValueFunct3D *BoundValue_Intl, TAuxParam3D *aux,
                            BoundCondFunct2D *SurfBDCond, BoundValueFunct2D * SurfBDValue )
 {

  char Name[] = "IntenalSystem";
  char Description[] = "InternalSystem";
  char Cstring[] = "U";
  char PosString[] = "PhysicalSpaceCoord";

  #ifdef _MPI
     if(SOLVER == DIRECT)
     {
      SQMATRICES_RTE[0] = sqmatrixM[N_Levels-1];
      TDS = new TParDirectSolver(ParComm[N_Levels-1],NULL,SQMATRICES_RTE,NULL);
    }
 #endif

 #ifdef _OMPONLY
    if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
    {
      DS = new TParDirectSolver(sqmatrixM[N_Levels-1]);
    }
 #endif 

  int i, N_SqMatrices, N_Rhs; 
  
   BoundaryConditions[0] = BoundCond;
   BoundaryConditions[1] = BoundCond;
   BoundaryConditions[2] = BoundCond;
   BoundaryConditions[3] = BoundCond;   
   BoundaryConditions[4] = BoundCond_Intl;
   
   BoundaryValues[0] = BoundValue;
   BoundaryValues[1] = BoundValue;
   BoundaryValues[2] = BoundValue;
   BoundaryValues[3] = BoundValue;         
   BoundaryValues[4] = BoundValue_Intl;

   TDiscreteForm3D *DiscreteFormBRhs_Galerkin, *DiscreteFormRhs_Galerkin, *DiscreteFormBRhs_SUPG, *DiscreteFormRhs_SUPG; 
    
   InitializeDiscreteFormRTE(DiscreteFormBRhs_Galerkin, DiscreteFormRhs_Galerkin, DiscreteFormBRhs_SUPG,
                             DiscreteFormRhs_SUPG, BilinearCoeffs);

     switch(Disctype)
      {
       case GALERKIN:
       case LOCAL_PROJECTION:
            DiscreteFormARhs = DiscreteFormBRhs_Galerkin;
            DiscreteRhs = DiscreteFormRhs_Galerkin;
       break;
       case SUPG:
      //  case SDFEM:    
            DiscreteFormARhs = DiscreteFormBRhs_SUPG;
            DiscreteRhs = DiscreteFormRhs_SUPG;
       break;

      
       default:
             OutPut("Unknown or not yet implemented DISCTYPE" << endl);
             exit(4711);;
      }  

    // initialize the assemble 
    if(aux==NULL)
     { aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); }
    
     for(i=Start_Level;i<N_Levels;i++)
     {  
      fesp[0] = FeSpaces[i];
      ferhs[0] = FeSpaces[i];  
      ferhs[1] = FeSpaces[i];  
      ferhs[2] = FeSpaces[i];  
      ferhs[3] = FeSpaces[i];  

      RHSs[0] = RhsArray[i];
      RHSs[1] = Rhs1Array[i];
      RHSs[2] = Rhs2Array[i];
      RHSs[3] = Rhs3Array[i];

      // A and B matrices Matrix     
      SQMATRICES_RTE[0] = sqmatrixM[i];
      SQMATRICES_RTE[1] = sqmatrixA[i];      
      SQMATRICES_RTE[2] = B1[i];
      SQMATRICES_RTE[3] = B2[i];
      SQMATRICES_RTE[4] = B3[i]; 
      N_SqMatrices = 5;
      N_Rhs = 4;

      if (Disctype==SUPG)  
      {
        SQMATRICES_RTE[5] = S1_SUPG[i];
        SQMATRICES_RTE[6] = S2_SUPG[i];
        SQMATRICES_RTE[7] = S3_SUPG[i];
        SQMATRICES_RTE[8] = S1S2_SUPG[i];
        SQMATRICES_RTE[9] = S1S3_SUPG[i];
        SQMATRICES_RTE[10] = S2S3_SUPG[i];   
        SQMATRICES_RTE[11] = M1_SUPG[i];   
        SQMATRICES_RTE[12] = M2_SUPG[i];   
        SQMATRICES_RTE[13] = M3_SUPG[i];   

        N_SqMatrices = 14;

        for(int j=0;j<9;j++)
         {
          RHSs[4+j] = Rhs_SUPG[j][i];   
          ferhs[4+j] = FeSpaces[i];  
          BoundaryConditions[4+j] = BoundCond;   
          BoundaryValues[4+j] = BoundValue;
         }

        N_Rhs = 13;
      }
 
      AMatRhsAssemble[i] = new TAssembleMat3D(1, fesp, N_SqMatrices, SQMATRICES_RTE, 0, NULL, N_Rhs, RHSs, ferhs, 
                               DiscreteFormARhs, BoundaryConditions, BoundaryValues, aux);
      AMatRhsAssemble[i]->Init(); 

      RhsAssemble[i] =  new TAssembleMat3D(1, fesp, 0, NULL, 0, NULL, N_Rhs, RHSs, ferhs, 
                               DiscreteRhs, BoundaryConditions, BoundaryValues, aux);
      RhsAssemble[i]->Init(); 
      //setup the multigrid solver
      if(SOLVER==GMG)
       {
 #ifdef _MPI  
        MGLevel = new TMGLevel3D(i, SQMATRICES_RTE[0], RHSs[0], SolArray[i], ParComm[i], ParMapper[i], N_aux, NULL);
 #else
        MGLevel = new TMGLevel3D(i, SQMATRICES_RTE[0], RHSs[0], SolArray[i], N_aux, NULL);
 #endif
        MG->AddLevel(MGLevel);
       }  
     } // for(i=Star 

   /** by default the intenal space is solved at the "N_Levels"th level of the physical space */
   SpatialSpace_Level = N_Levels-1;

   N_SpatialPts = FeSpaces[SpatialSpace_Level]->GetN_DegreesOfFreedom();
   SpatialPts = new double[3*N_SpatialPts];
   SpatialPos = new TFEVectFunct3D(FeSpaces[SpatialSpace_Level], PosString, PosString, SpatialPts, N_SpatialPts, 3);

   SpatialPos->GridToData();
   Spatial_X = SpatialPts;
   Spatial_Y = SpatialPts+N_SpatialPts;
   Spatial_Z = SpatialPts+N_SpatialPts;           

  //  cout << " N_SpatialPts : " << N_SpatialPts <<endl;
  //=========================================================================
  /** Init for internal spaces */
  //=========================================================================  

  /** set data for multigrid */
  ORDER_INTL = TDatabase::ParamDB->ANSATZ_ORDER_INTL;

//=========================================================================
// construct all internal finite element spaces
//=========================================================================
  // fe (bulk) space for internal equation, needed for interpolation and output     
  FeSpaces_Intl3D =  new TFESpace3D(Coll_Intl, Name, Description, BoundCond_Intl, Non_USpace, ORDER_INTL+1);  
  N_DoF_Intl3D = FeSpaces_Intl3D->GetN_DegreesOfFreedom();
  Sol_Intl3D = new double[N_DoF_Intl3D];

  FeFunction_Intl3D  = new TFEFunction3D(FeSpaces_Intl3D, Cstring, Cstring, Sol_Intl3D, N_DoF_Intl3D);  

  std::ostringstream os;
  os << " ";
  Output = new TOutput3D(1, 1, 1, 1, Domain_Intl);
  os.seekp(std::ios::beg);
  Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
  Output->AddFEFunction(FeFunction_Intl3D);  

  os.seekp(std::ios::beg);
  os << TDatabase::ParamDB->VTKBASENAME<<".vtk" << ends;
  Output->WriteVtk(os.str().c_str());

 // fe (surface) space for internal equation
  /** physcal space is solved on the surface of the internal domain */
  // cout << "Number of Surface Cells: " << N_SurfCells<<endl;  

  FeSpaces_Intl = new TFESpace2D(Surf_Coll, Name, Description, SurfBDCond, ORDER_INTL, NULL);
  N_DoF_Intl = FeSpaces_Intl->GetN_DegreesOfFreedom();
  cout << "N_Surf_Dof : " << N_DoF_Intl <<endl;
  cout << "ORDER_INTL : " << ORDER_INTL <<endl;
 
  sol_intl = new double[N_DoF_Intl*N_SpatialPts];
  rhs_intl = new double[N_DoF_Intl*N_SpatialPts];
  BE_MatScale = new double[N_DoF_Intl];
  //get the coordinates of the internal points (surface points), needed for physical space convection term
  InternalPts = new double[3*N_DoF_Intl];
  
  FeSpaces_Intl3D->GetSurfDOFPosition(N_SurfCells, Cell_IDX, Joint_IDX, InternalPts);

  //  cout << " N_DoF_Intl : " << N_DoF_Intl <<endl;
  //  for(i=0;i<N_DoF_Intl;i++)
  //  cout << " InternalPts : " << InternalPts[3*i] << " " << InternalPts[3*i+1]<< " " << InternalPts[3*i+2]<<endl;

  // cout << "SystemTCD3D_3C.C internal Init done! " <<endl;
 
} // Init


void TSystemRTE::Interpolate(DoubleFunct3D *Exact, double *Sol)
{
 double *sol = SpatialFEFunctions[SpatialSpace_Level]->GetValues();

  for(int i=0; i<N_DoF_Intl; ++i)
   {
    TDatabase::ParamDB->P11 = InternalPts[3*i];
    TDatabase::ParamDB->P12 = InternalPts[3*i+1];
    TDatabase::ParamDB->P13 = InternalPts[3*i+2];
    SpatialFEFunctions[SpatialSpace_Level]->Interpolate(Exact);
      
    memcpy(Sol+(i*N_SpatialPts), sol, N_SpatialPts*SizeOfDouble); 
   } 
}
 
void TSystemRTE::AssembleABRhs()
 {
  //this is set to true for direct solver factorization
  factorize = true;
  
  int i, N_DOF_low, N_Active;

  for(i=Start_Level;i<N_Levels;i++)
  {    
   /** reset the matrix and rhs */
   AMatRhsAssemble[i]->Reset(); 

   // assemble
   AMatRhsAssemble[i]->Assemble3D();
                     
  }//   for(i=Start_Level;i<N_Le  
    
} // TSystemMatScalar3D::AssembleARhs 


void TSystemRTE::AssembleRhs()
 {
  int i, N_DOF_low, N_Active;

  memcpy(rhs1old, Rhs1Array[N_Levels-1], N_DOF*SizeOfDouble);  
  memcpy(rhs2old, Rhs2Array[N_Levels-1], N_DOF*SizeOfDouble);  
  memcpy(rhs3old, Rhs3Array[N_Levels-1], N_DOF*SizeOfDouble);  
  

  for(i=Start_Level;i<N_Levels;i++)
  {    
   /** reset the matrix and rhs */
   RhsAssemble[i]->Reset(); 

   // assemble 
   RhsAssemble[i]->Assemble3D();                  
  }//   for(i=Start_Level;i<N_Le  


 } // AssembleRhs()

void TSystemRTE::AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol, 
                                 double b1, double b2, double b3, CoeffFct3D *SystemRhs
#ifdef _MPI
                                 , double **Rhs_array
#endif
                                              )
 {
     int i, N_Active;
     double tau, param[3];
    // double diffusion = TDatabase::ParamDB->PE_NR;

     if(SystMatAssembled)
      {
       OutPut("System has to be restored before calling AssembleSystMat! " <<endl);
       exit(0);
      }

     N_Active =  FeSpaces[N_Levels-1]->GetActiveBound();     
     tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;  
    
     memset(B, 0, N_DOF*SizeOfDouble); 

     /** old rhs multiplied with current subtime step and theta3 on B */
     if(fabs(TDatabase::TimeDB->THETA3)>0)
      {
       Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3,  oldrhs, B);    
       Daxpy(N_Active, b1*tau*TDatabase::TimeDB->THETA3,  rhs1old, B);   
       Daxpy(N_Active, b2*tau*TDatabase::TimeDB->THETA3,  rhs2old, B);   
       Daxpy(N_Active, b3*tau*TDatabase::TimeDB->THETA3,  rhs3old, B);    

      cout<<"Old Rhs..Array in AssembleRhs() need to be stored befrore using this time-disc! " <<endl;
       exit(0);

      }

     /** add rhs from current sub time step to rhs array B */
     if(fabs(TDatabase::TimeDB->THETA4)>0)
     {
      RHSs[0] = rhs;
      RHSs[1] = Rhs1Array[N_Levels-1];
      RHSs[2] = Rhs2Array[N_Levels-1];
      RHSs[3] = Rhs3Array[N_Levels-1];

      if (Disctype==SUPG)
      {
        for(int j=0;j<9;j++)
          RHSs[4+j] = Rhs_SUPG[j][N_Levels-1];
      }

      param[0] = b1;
      param[1] = b2;      
      param[2] = b3;      
      SystemRhs(N_Active, param, B, NULL, RHSs, NULL);
     }

      /** M + \delta_K (U, S\cdot\grad V) */
     if (Disctype==SUPG)
      {
       //time consistent SUPG
       MatAdd(sqmatrixM[N_Levels-1], M1_SUPG[N_Levels-1], b1);
       MatAdd(sqmatrixM[N_Levels-1], M2_SUPG[N_Levels-1], b2);
       MatAdd(sqmatrixM[N_Levels-1], M3_SUPG[N_Levels-1], b3);
      }

     /** M = M + (- tau*THETA2)A */
     if(fabs(TDatabase::TimeDB->THETA2)>0)
     {
      MatAdd(sqmatrixM[N_Levels-1], sqmatrixA[N_Levels-1], -tau*TDatabase::TimeDB->THETA2);
      MatAdd(sqmatrixM[N_Levels-1], B1[N_Levels-1], - b1*tau*TDatabase::TimeDB->THETA2);
      MatAdd(sqmatrixM[N_Levels-1], B2[N_Levels-1], - b2*tau*TDatabase::TimeDB->THETA2);
      MatAdd(sqmatrixM[N_Levels-1], B3[N_Levels-1], - b3*tau*TDatabase::TimeDB->THETA2);

      gammab0 = -tau*TDatabase::TimeDB->THETA2;  // set current factor of steady statmatrix
      gammab1 = -b1*tau*TDatabase::TimeDB->THETA2;  // set current factor of steady statmatrix
      gammab2 = -b2*tau*TDatabase::TimeDB->THETA2;  // set current factor of steady state matrix
      gammab3 = -b3*tau*TDatabase::TimeDB->THETA2;  // set current factor of steady state matrix

      if (Disctype==SUPG)
      {
       MatAdd(sqmatrixM[N_Levels-1], S1_SUPG[N_Levels-1], -b1*b1*tau*TDatabase::TimeDB->THETA2);
       MatAdd(sqmatrixM[N_Levels-1], S2_SUPG[N_Levels-1], -b2*b2*tau*TDatabase::TimeDB->THETA2);
       MatAdd(sqmatrixM[N_Levels-1], S3_SUPG[N_Levels-1], -b3*b3*tau*TDatabase::TimeDB->THETA2);
 
       MatAdd(sqmatrixM[N_Levels-1], S1S2_SUPG[N_Levels-1], - b1*b2*tau*TDatabase::TimeDB->THETA2);
       MatAdd(sqmatrixM[N_Levels-1], S1S3_SUPG[N_Levels-1], - b1*b3*tau*TDatabase::TimeDB->THETA2);
       MatAdd(sqmatrixM[N_Levels-1], S2S3_SUPG[N_Levels-1], - b2*b3*tau*TDatabase::TimeDB->THETA2); 

       gammab4 = -b1*b1*tau*TDatabase::TimeDB->THETA2;  // set current factor of steady statmatrix
       gammab5 = -b2*b2*tau*TDatabase::TimeDB->THETA2;  // set current factor of steady state matrix
       gammab6 = -b3*b3*tau*TDatabase::TimeDB->THETA2;  // set current factor of steady state matrix
       gammab7 = -b1*b2*tau*TDatabase::TimeDB->THETA2;  // set current factor of steady statmatrix
       gammab8 = -b1*b3*tau*TDatabase::TimeDB->THETA2;  // set current factor of steady state matrix
       gammab9 = -b2*b3*tau*TDatabase::TimeDB->THETA2;  // set current factor of steady state matrix

      }


     }
     else
     {
      gammab0 = 0.;     
      gammab1 = 0.;   
      gammab2 = 0.;   
      gammab3 = 0.;   
      gammab4 = 0.;
      gammab5 = 0.; 
      gammab6 = 0.; 
      gammab7 = 0.; 
      gammab8 = 0.; 
      gammab9 = 0.;  
     }
    /** defect = M * oldsol */
     memset(defect, 0, N_DOF*SizeOfDouble);  
     MatVectActive(sqmatrixM[N_Levels-1], oldsol, defect); 
     //cout << "defect " << Ddot(N_Active, sol, sol)<< endl;      
    
    /** B:= B + defec  */
     Daxpy(N_Active, 1, defect, B);

    /** set Dirichlet values */
     memcpy(B+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);  
     memcpy(sol+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
          
    //for rteSpatialINternal.h example only (hard coded)
    Dscal((N_DOF-N_Active), b1*b2*b3, B+N_Active);
    Dscal((N_DOF-N_Active), b1*b2*b3, sol+N_Active);     
 
    /** assemble the system matrix */
      for(i=Start_Level;i<N_Levels;i++)   
       {
        if(i==N_Levels-1)
         {
     
          MatAdd(sqmatrixM[i], sqmatrixA[i], -gammab0 + tau*TDatabase::TimeDB->THETA1);

          MatAdd(sqmatrixM[i], B1[i], -gammab1 + b1*tau*TDatabase::TimeDB->THETA1);
          MatAdd(sqmatrixM[i], B2[i], -gammab2 + b2*tau*TDatabase::TimeDB->THETA1);
          MatAdd(sqmatrixM[i], B3[i], -gammab3 + b3*tau*TDatabase::TimeDB->THETA1);    

          if (Disctype==SUPG)
          {
           MatAdd(sqmatrixM[i], S1_SUPG[i], -gammab4 + b1*b1*tau*TDatabase::TimeDB->THETA1);
           MatAdd(sqmatrixM[i], S2_SUPG[i], -gammab5 + b2*b2*tau*TDatabase::TimeDB->THETA1);
           MatAdd(sqmatrixM[i], S3_SUPG[i], -gammab6 + b3*b3*tau*TDatabase::TimeDB->THETA1);
 
           MatAdd(sqmatrixM[i], S1S2_SUPG[i], -gammab7 + b1*b2*tau*TDatabase::TimeDB->THETA1);
           MatAdd(sqmatrixM[i], S1S3_SUPG[i], -gammab8 + b1*b3*tau*TDatabase::TimeDB->THETA1);
           MatAdd(sqmatrixM[i], S2S3_SUPG[i], -gammab9 + b2*b3*tau*TDatabase::TimeDB->THETA1); 
          }
         }
        else
         { 
          MatAdd(sqmatrixM[i], sqmatrixA[i], tau*TDatabase::TimeDB->THETA1);

          MatAdd(sqmatrixM[i], B1[i],  b1*tau*TDatabase::TimeDB->THETA1);
          MatAdd(sqmatrixM[i], B2[i],  b2*tau*TDatabase::TimeDB->THETA1);
          MatAdd(sqmatrixM[i], B3[i],  b3*tau*TDatabase::TimeDB->THETA1); 

          if (Disctype==SUPG)
          {
           //time consistent SUPG on low MG levels
           MatAdd(sqmatrixM[i], M1_SUPG[i], b1);
           MatAdd(sqmatrixM[i], M2_SUPG[i], b2);
           MatAdd(sqmatrixM[i], M3_SUPG[i], b3);

           MatAdd(sqmatrixM[i], S1_SUPG[i], b1*b1*tau*TDatabase::TimeDB->THETA1);
           MatAdd(sqmatrixM[i], S2_SUPG[i], b2*b2*tau*TDatabase::TimeDB->THETA1);
           MatAdd(sqmatrixM[i], S3_SUPG[i], b3*b3*tau*TDatabase::TimeDB->THETA1);
 
           MatAdd(sqmatrixM[i], S1S2_SUPG[i], b1*b2*tau*TDatabase::TimeDB->THETA1);
           MatAdd(sqmatrixM[i], S1S3_SUPG[i], b1*b3*tau*TDatabase::TimeDB->THETA1);
           MatAdd(sqmatrixM[i], S2S3_SUPG[i], b2*b3*tau*TDatabase::TimeDB->THETA1); 
          }


        } 
         
 #ifdef _MPI  
        SQMATRICES_RTE[0] = sqmatrixM[i];  
 #endif
        } // for(i=Start_Leve

     gammab0 = tau*TDatabase::TimeDB->THETA1;
     gammab1 = b1*tau*TDatabase::TimeDB->THETA1;  // set current factor of steady state matrix
     gammab2 = b2*tau*TDatabase::TimeDB->THETA1;  // set current factor of steady state matrix
     gammab3 = b3*tau*TDatabase::TimeDB->THETA1;  // set current factor of steady state matrix

     if (Disctype==SUPG)
      {
       gammab4 = b1*b1*tau*TDatabase::TimeDB->THETA1; 
       gammab5 = b2*b2*tau*TDatabase::TimeDB->THETA1; 
       gammab6 = b3*b3*tau*TDatabase::TimeDB->THETA1; 
       gammab7 = b1*b2*tau*TDatabase::TimeDB->THETA1; 
       gammab8 = b1*b3*tau*TDatabase::TimeDB->THETA1; 
       gammab9 = b2*b3*tau*TDatabase::TimeDB->THETA1; 
       gammab10 = b1;
       gammab11 = b2;
       gammab12 = b3;
      }

// //have to shift this in pardirectsolver     
// #ifdef _OMPONLY     
//     if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
//       DS->AssembleMatrix(sqmatrixM[N_Levels-1]);
// #endif
     
  SystMatAssembled  = TRUE;
 // cout << "System is assembled "  << b1 <<endl;
} // AssembleSystMat

void TSystemRTE::RestoreMassMat()
 {
  int i;

   if(SystMatAssembled)
    {
      // restore the mass matrix
      for(i=Start_Level;i<N_Levels;i++)  
       {
        MatAdd(sqmatrixM[i], sqmatrixA[i], -gammab0);    

        MatAdd(sqmatrixM[i], B1[i], -gammab1);
        MatAdd(sqmatrixM[i], B2[i], -gammab2);
        MatAdd(sqmatrixM[i], B3[i], -gammab3);    

        if (Disctype==SUPG)
          {

           //time consistent SUPG
           MatAdd(sqmatrixM[i], M1_SUPG[i], -gammab10);
           MatAdd(sqmatrixM[i], M2_SUPG[i], -gammab11);
           MatAdd(sqmatrixM[i], M3_SUPG[i], -gammab12);
  
           MatAdd(sqmatrixM[i], S1_SUPG[i], -gammab4 );
           MatAdd(sqmatrixM[i], S2_SUPG[i], -gammab5 );
           MatAdd(sqmatrixM[i], S3_SUPG[i], -gammab6 );
 
           MatAdd(sqmatrixM[i], S1S2_SUPG[i], -gammab7 );
           MatAdd(sqmatrixM[i], S1S3_SUPG[i], -gammab8 );
           MatAdd(sqmatrixM[i], S2S3_SUPG[i], -gammab9 ); 
        }        

       }
     
      gammab0 = 0;
      gammab1 = 0.;
      gammab2 = 0.;
      gammab3 = 0.;
      if (Disctype==SUPG)
       {
        gammab4 = 0.;
        gammab5 = 0.; 
        gammab6 = 0.; 
        gammab7 = 0.; 
        gammab8 = 0.; 
        gammab9 = 0.;  
        gammab10 = 0;
        gammab11 = 0;
        gammab12 = 0;
       }

    SystMatAssembled  = FALSE;
   }
  else
   {
     cout << "System is not assembled to restore " <<endl;
     exit(0);
   }

 }

void TSystemRTE::Solve(double *sol)
 {  
  switch(SOLVER)
   {
//       case AMG_SOLVE:
//          Solver(sqmatrixM[N_Levels-1], B, sol);
//       break;

//       case GMG:
//         if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
//          {
//           memcpy(Itmethod_sol, sol, N_DOF*SizeOfDouble);
//           memcpy(Itmethod_rhs, B, N_DOF*SizeOfDouble);
//          }
//         else
//          {
//           Itmethod_sol = sol;
//           Itmethod_rhs = B;
//          }
      
//          /** solve linear system */
//         Itmethod->Iterate(sqmatrices, NULL, Itmethod_sol, Itmethod_rhs);
// #ifdef _MPI
//     if(TDatabase::ParamDB->SC_SMOOTHER_SCALAR==6)
//          ParComm[N_Levels-1]->CommUpdateH2(Itmethod_sol);
// #endif
//         if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
//          {
//           memcpy(sol, Itmethod_sol, N_DOF*SizeOfDouble);
//          }
//       break;

       case DIRECT:
 #ifdef _MPI
 	TDS->Solve(sol, B, factorize);
 // 	exit(0);
 #endif

 #ifdef _OMPONLY
 	if(TDatabase::ParamDB->DSType == 1)
 	  DS->Solve(sol, B, factorize);
 	else{
 	  OutPut("Select Proper Solver" << endl);
 	  exit(0);
 	}
 #endif

 #ifdef _SEQ
    DirectSolver(sqmatrixM[N_Levels-1], B, sol);

   // cout << "seq direct solver " <<endl;
 #endif
//this is set to false for direct solver factorization
    factorize = false;
    break;      
 
  default:
    OutPut("Unknown Solver" << endl);
    exit(4711);;
   }
     
 } //solve

void TSystemRTE::GetQuadIntegral(int i, int &N_Points, double * &weights, double *x, double *y, double *z,
                                 double *n1, double *n2, double *n3)
{
 int k, Cell_No, Joint_No;
 double *p1, *p2, xi, eta, zeta;
 double a1, a2, a3, b1, b2, b3;

 TBaseCell *Cell;
 TJoint *Joint;
 bool Isoparametric;
 RefTrans3D RefTrans;
 TRefTrans3D *F_K;
 QuadFormula2D QuadFormula;
 TQuadFormula2D *qf2d;

   Cell_No = Cell_IDX[i];
   Joint_No = Joint_IDX[i];
   Cell = Coll_Intl->GetCell(Cell_No);
   Joint = Cell->GetJoint(Joint_No);

   if (Joint->GetType() == IsoBoundFace)
    { Isoparametric = TRUE; }
    
 	  RefTrans = TetraAffin;
	  F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
	  ((TTetraAffin*) F_K)->SetCell(Cell);
	  QuadFormula = Gauss3Tria;
    qf2d = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
    qf2d->GetFormulaData(N_Points, weights, p1, p2);

    // compute normals for all quad points 
    for (k=0;k<N_Points;++k)
    {
	    switch (Joint_No)
	    {
	      case 0:
		      xi = p1[k];
		      eta = p2[k];
		      zeta = 0;
		    break;
	      case 1:
		      xi = p2[k];
		      eta = 0;
		      zeta = p1[k];
		    break;
	      case 2:
		      xi = p1[k];
		      eta = 1-p1[k]-p2[k];
		      zeta = p2[k];
		    break;
	      case 3:
		      xi = 0;
		      eta = p1[k];
		      zeta = p2[k];
		    break;
	    }	 

	   ((TTetraAffin*) F_K)->GetTangentVectors(Joint_No, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
      n1[k] = (a2*b3 - a3*b2);
      n2[k] = (a3*b1 - a1*b3);
      n3[k] = (a1*b2 - a2*b1);

      ((TTetraAffin*) F_K)->GetOrigFromRef(xi, eta, zeta, x[k], y[k], z[k]);      
    } //  for (k=0;k<N_Points;++



} // void TSystemRTE::GetQuadIntegral(i


void TSystemRTE::AssembleIntlMat(int CellNo, double tau, double *XPos, double &fact, double *InternalScaling)
{
 int i, k, kk, l, m, Cell_No, Joint_No, N_Points;
 double *weights, *weightsk, p1, *p2, xi, eta, zeta;
 double a1, a2, a3, b1, b2, b3, sigma_a, sigma_s;
 double val, len, len1, surf_area=0., surf_areak, stiff;
 double n1[MaxN_QuadPoints_2D], n2[MaxN_QuadPoints_2D], n3[MaxN_QuadPoints_2D];
 double x[MaxN_QuadPoints_2D], y[MaxN_QuadPoints_2D], z[MaxN_QuadPoints_2D];
 double nk1[MaxN_QuadPoints_2D], nk2[MaxN_QuadPoints_2D], nk3[MaxN_QuadPoints_2D];
 double xk[MaxN_QuadPoints_2D], yk[MaxN_QuadPoints_2D], zk[MaxN_QuadPoints_2D];

 double Coords[9], InnVal[1], mass;

 TBaseCell *Cell;
 TJoint *Joint;
 bool Isoparametric;
 RefTrans3D RefTrans;
 TRefTrans3D *F_K;
 QuadFormula2D QuadFormula;
 TQuadFormula2D *qf2d;

 sigma_a = TDatabase::ParamDB->P1;
 sigma_s = TDatabase::ParamDB->P2;
 
//  memset(BE_MatScale, 0, N_DoF_Intl*SizeOfDouble); 
 memset(InternalScaling, 0, N_SpatialPts*SizeOfDouble); 
 

 if(N_DoF_Intl!=N_SurfCells)
  {
    cout<<"Only Po is allowed for surface mesh"<<endl;
    exit(0);
  }


 i=CellNo;
//  for(i=0;i<N_SurfCells; i++)
 {
    this->GetQuadIntegral(i, N_Points, weights, x, y, z, n1, n2, n3);

    mass = 0.;
    stiff = 0.;
    for (k=0;k<N_Points;++k)
    {
      len = sqrt(n1[k]*n1[k] + n2[k]*n2[k] + n3[k]*n3[k]);
      len *= weights[k]; /** X * n1*/;

      mass  += len;
      stiff += sigma_a*len;      

      Coords[0] = x[k];
      Coords[1] = y[k];      
      Coords[2] = z[k];
      //compute rhs 
      for(m=0;m<N_SpatialPts; ++m)
      {
       Coords[6] = XPos[m]; //spatial x
       Coords[7] = XPos[N_SpatialPts+m];      
       Coords[8] = XPos[2*N_SpatialPts+m];   

      //  get rhs     
       GetRhs(6, Coords, InnVal);
       
       InternalScaling[m] +=len*InnVal[0]; //  rhs +=len*InnVal[0];
      }
      //  since \int_S^2(kernel) = 1 and in dG(0), nalways 1, so no need to compute
      val = 1.0;
      // val = 0;
      // // integral (double) of the kernal function
      // for(kk=0;kk<N_SurfCells; kk++)
      //  {
      //   this->GetQuadIntegral(kk, N_Points, weightsk, xk, yk, zk, nk1, nk2, nk3);

      //   for (l=0;l<N_Points;++l)
      //    {
      //     Coords[3] = xk[l];
      //     Coords[4] = yk[l];      
      //     Coords[5] = zk[l];

      //     len1 = weightsk[l]*sqrt(nk1[l]*nk1[l] + nk2[l]*nk2[l] + nk3[l]*nk3[l]);

      //     //get kernel values 
      //     GetKernel(6, Coords, InnVal);

      //     //  stiff += val*len1*InnVal[0];
      //      val += len1*InnVal[0];
      //     } //  for (int l=0;l<N_Points
      //  } // 

       
      // OutPut(" surf_areak: " <<  val <<endl);

      stiff -=sigma_s*len*val;  
      // surf_area += len;
    }//for (int k=

   // rhs scaling factor in backwar Euler
  // BE_MatScale[i] = (1. - tau*TDatabase::TimeDB->THETA2*stiff)/ (1. + tau*TDatabase::TimeDB->THETA1*stiff);

   fact = 1./(1. +tau*TDatabase::TimeDB->THETA1*stiff);
  //only backward Euler
  for(m=0;m<N_SpatialPts; ++m)
   {
     InternalScaling[m] *= tau*fact;
    //  OutPut(fact << " Mat: " <<  InternalScaling[m] <<endl);     
   }

 } // for(i=0;i
// tau*TDatabase::TimeDB->THETA1
//  OutPut(4.*Pi<< " Surf Area: " << surf_area <<endl);
//  exit(0);

} // void TSystemRTE::AssembleIntlM


void TSystemRTE::GetErrors(DoubleFunctVect *Exact, TFEFunction3D *ScalarFunction, double *sol_all, double &l2)
{
 int i, j, k, l, m, P, N_Cells, N_LDof, N_LocalUsedElements, *N_BaseFunct, N_U;
 int *BeginIndex, *GlobalNumbers, *BeginIndex_Intl, *GlobalNumbers_Intl;
 int N_Points, N_LocalDOFs, *DOF, N_BaseFunct_Intl, *DOF_Intl, N_LinePoints, N_Sets=1;
  
 double *weights, *xi, *eta, *zeta;
 double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
 double AbsDetjk[MaxN_QuadPoints_3D];
 double *Sol_QuadIntl, *sol;
 double **OrigFEValues, *Orig, value, InternalError;
 double S1[MaxN_QuadPoints_2D], S2[MaxN_QuadPoints_2D], S3[MaxN_QuadPoints_2D];
 double **origvaluesD0, **origvaluesD1, *orgD0, *orgD1, Mult;
 double  Pos[6], Exactsol[6], *solIntl, u;
 double errors = 0.;
 
 TBaseCell *cell, *Cell_Intl;
 TCollection *Coll;
 TFESpace3D *FESpace;  
 BaseFunct3D BaseFunct, *BaseFuncts;
 TFE3D *Element;
 FE3D LocalUsedElements[1], CurrentElement;
 FE1D FEId_Intl;
 TFE1D *Element_Intl;
 TBaseFunct1D *bf_Intl;
 BaseFunct1D BaseFunct_ID_Intl, BaseFunct_Intl[1];
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans3D *F_K;
 RefTrans3D RefTrans;
 double **uref;

 int Cell_No, Joint_No, lIntl, ii, N_PointsIntl;
 TBaseCell *CellIntl3D;
 TJoint *Joint;
 FE3D FEIdIntl;
 BF3DRefElements RefElementIntl;
 QuadFormula2D QuadFormula;
 TQuadFormula2D *qf;
 double *Weights, *p1, *p2, *WeightsIntl;
 double a1, a2, a3, b1, b2, b3, n1, n2, n3, len;
 double surfArea, volume;
 bool *SecondDer;
 bool Needs2ndDer[1];

  FESpace = ScalarFunction->GetFESpace3D();
  N_U = ScalarFunction->GetLength();  
  BeginIndex = FESpace->GetBeginIndex();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex_Intl = FeSpaces_Intl->GetBeginIndex();
  GlobalNumbers_Intl = FeSpaces_Intl->GetGlobalNumbers();
  N_LDof = FeSpaces_Intl->GetN_DegreesOfFreedom();
  N_LocalUsedElements = 1;
  SecondDer = new bool[1];
  SecondDer[0] = FALSE;
  Needs2ndDer[0] = FALSE;

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
 
  //first check how many quad pts  in the cell
  cell = Coll->GetCell(0);
  LocalUsedElements[0] = FESpace->GetFE3D(0, cell);
  Element = TFEDatabase3D::GetFE3D(LocalUsedElements[0]);
  TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, Coll, cell, (bool *)SecondDer,
                         N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
  
  Sol_QuadIntl = new double[N_Points*N_DoF_Intl];
  
  // volume = 0;

  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    LocalUsedElements[0] = FESpace->GetFE3D(i, cell);
    Element = TFEDatabase3D::GetFE3D(LocalUsedElements[0]);

    // ====================================================================
    // calculate values on original element
    // ====================================================================
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, Coll, cell, (bool *)SecondDer,
                           N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
 
    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace->GetFE3D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_LocalDOFs = N_BaseFunct[CurrentElement];
    DOF = GlobalNumbers+BeginIndex[i];

    OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);

    // find values for all internal levels at all quad points in ith cell 
    memset(Sol_QuadIntl, 0, N_DoF_Intl*N_Points*SizeOfDouble);
    
    for(j=0; j<N_DoF_Intl; j++)
    {
      sol = sol_all+j*N_U;

      for(k=0; k<N_Points; k++)
      {
        Orig = OrigFEValues[k];
        value = 0.;

        for(l=0; l<N_LocalDOFs; l++)
          value += sol[DOF[l]]*Orig[l];

        Sol_QuadIntl[k*N_DoF_Intl  +  j] = value;
      }
    }      //for(j=0; j<N_DoF_Intl; j++)    

 
    if(N_SurfCells != N_DoF_Intl)
     {
      cout << "Only Po is implemented for Internal" <<endl;
      cout<< N_SurfCells << " N_DoF_Intl " << N_DoF_Intl << endl;
      exit(0);
      /** can be extended for P_k^{disc} elements as well */
     } 

    for(j=0; j<N_Points; j++)
     {
      InternalError = 0;
      Pos[0] = X[j];
      Pos[1] = Y[j];
      Pos[2] = Z[j];
      // surfArea = 0.;
      solIntl = Sol_QuadIntl+(j*N_SurfCells);
      
      for(ii=0;ii<N_SurfCells; ++ii)
       {
        /** since the internal is P0 space, so sol is constant on each surf cell */
        u = solIntl[ii];
   
        Cell_No = Cell_IDX[ii];
        Joint_No = Joint_IDX[ii];
        CellIntl3D = Coll_Intl->GetCell(Cell_No);
        Joint = CellIntl3D->GetJoint(Joint_No);
        FEIdIntl = FeSpaces_Intl3D->GetFE3D(Cell_No, CellIntl3D);
        RefElementIntl = TFEDatabase3D::GetRefElementFromFE3D(FEIdIntl);
        lIntl = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEIdIntl);

        switch(RefElementIntl)
         {
          case BFUnitTetrahedron :
          RefTrans = TetraIsoparametric;
	        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
	        ((TTetraIsoparametric*) F_K)->SetApproximationOrder(1);
	        ((TTetraIsoparametric*) F_K)->SetCell(CellIntl3D);
          break;

          case BFUnitHexahedron :
            Error("Nothing is implemented for Hexahedron!" << endl);
            exit(-1);
          break;
        } // endswitch

        QuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*lIntl);
        qf = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
        qf->GetFormulaData(N_PointsIntl, WeightsIntl, p1, p2);

        switch(RefElementIntl)
         {
          case BFUnitTetrahedron :
            ((TTetraIsoparametric *)F_K)->GetOrigBoundFromRef(Joint_No, N_PointsIntl, p1, p2, S1, S2, S3);
          break;

          case BFUnitHexahedron:
           Error("Nothing is implemented for Hexahedron!" << endl);
           exit(-1);
          break;
         } // endswitch

        for(k=0;k<N_PointsIntl;k++)
         {
          Pos[3] = S1[k];
          Pos[4] = S2[k];
          Pos[5] = S3[k];

          switch(RefElementIntl)
           {
            case BFUnitTetrahedron :
             ((TTetraIsoparametric *)F_K)->GetTangentVectors(Joint_No, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
            break;

            case BFUnitHexahedron:
            Error("Nothing is implemented for tetrahedra!" << endl);
            exit(-1);
           break;
          } // endswitch

          n1 = a2*b3 - a3*b2;
          n2 = a3*b1 - a1*b3;
          n3 = a1*b2 - a2*b1;

          len = sqrt(n1*n1 + n2*n2 + n3*n3);
          Mult = len*WeightsIntl[k];

          Exact(Pos, Exactsol);
          InternalError += Mult*(Exactsol[0]-u)*(Exactsol[0]-u);
          // InternalError += Mult;
          // surfArea +=  Mult;
         }//  for(k=0;k<N_PointsIntl;k++)
       } //   for(ii=0;ii<N_SurfCells; ++ii''

      // compute the error for all spatial quad points  
      Mult = weights[j]*AbsDetjk[j];
      errors += Mult*InternalError;  
    } // for(j=0; j<N_Points; j++)

  } // for(i=0; i<N_Cells; i++)

#ifdef _SMPI
  double recvbuf;

  MPI_Reduce(&errors, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, TDatabase::ParamDB->Comm );
  l2 = sqrt(recvbuf);
#else
  l2 = sqrt(errors);  
#endif



  //  OutPut( "L2: " << errors  << endl);
  //  cout <<   "Volume : " << volume <<endl;

  delete [] Sol_QuadIntl;
  // delete [] Sol_AllL;
  // delete [] InternalError;
   
} // void TSystemRTE::GetErrors

// double TSystemTCD3D::GetResidual(double *sol)
// {
//   double residual_scalar=0.0;
  
//   if(SystMatAssembled)
//    {
//     memset(defect, 0, N_DOF*SizeOfDouble);         
//     ScalarDefect(sqmatrixM[N_Levels-1], sol, B, defect, residual_scalar);
    
// #ifdef _MPI 
//     residual_scalar = 0.0;
//     double sum =0.;
//     int i,rank;
//     MPI_Comm_rank(Comm, &rank);
//     int *master = ParComm[N_Levels-1]->GetMaster();
//     for(i=0;i<N_DOF;i++)
//     {
//       if(master[i]!=rank)    continue;
//       residual_scalar += defect[i]*defect[i];
//     }
//    MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
//    residual_scalar = sqrt(sum);
// #endif
//    }
//   else
//    {
//     OutPut("Assemble the System Matrix before calculating the GetResidual" << endl);
//     exit(4711);;   
//    }
//    return residual_scalar;    
// }



