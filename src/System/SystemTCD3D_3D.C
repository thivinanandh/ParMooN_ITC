/** ************************************************************************ 
* @brief     source file for TSystemTCD3D_3D
* @author    Sashikumaar Ganesan
* @date      06.10.2019
* @History 
 ************************************************************************  */
 #include <Database.h>
 #include <SystemTCD3D_3D.h>
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

#include <stdlib.h>
#include <string.h>

TSystemTCD3D_3D::TSystemTCD3D_3D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver, int intspacedim)
                       :TSystemTCD3D(N_levels, fespaces, sol, rhs, disctype, solver)
{
 char  *GEO, *PRM;

 InternalSpaceDim = intspacedim;
 //=========================================================================
 // Internal space
 // construct FESpace and matrices for internal space matrices
 //=========================================================================
  Domain_Intl = new TDomain();  
  PRM = TDatabase::ParamDB->BNDFILE_INTL;
  GEO = TDatabase::ParamDB->GEOFILE_INTL;

  //generate mesh
  Domain_Intl->Init(PRM, GEO); // ParMooN  build-in Geo mesh

  LEVELS_INTL = TDatabase::ParamDB->LEVELS_INTL;
  if(TDatabase::ParamDB->SOLVER_TYPE_INTL==DIRECT)
  {
    TDatabase::ParamDB->UNIFORM_STEPS_INTL += (LEVELS_INTL-1);
    LEVELS_INTL = 1;
  }

  for(int i=0;i<TDatabase::ParamDB->UNIFORM_STEPS_INTL;i++)
    Domain_Intl->RegRefineAll();

  Coll_Intl=Domain_Intl->GetCollection(It_Finest, 0);


} // constructor


TSystemTCD3D_3D::~TSystemTCD3D_3D()
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


 void TSystemTCD3D_3D::Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue,
                            BoundCondFunct3D *BoundCond_Intl, BoundValueFunct3D *BoundValue_Intl, TAuxParam3D *aux )
 {

  char Name[] = "IntenalSystem";
  char Description[] = "InternalSystem";
  char Cstring[] = "U";
  char PosString[] = "PhysicalSpaceCoord";

  #ifdef _MPI
     if(SOLVER == DIRECT)
     {
      SQMATRICES[0] = sqmatrixM[N_Levels-1];
      TDS = new TParDirectSolver(ParComm[N_Levels-1],NULL,SQMATRICES,NULL);
    }
 #endif

 #ifdef _OMPONLY
    if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
    {
      DS = new TParDirectSolver(sqmatrixM[N_Levels-1]);
    }
 #endif 

  int i; 
  
   BoundaryConditions[0] = BoundCond;
   BoundaryConditions[1] = BoundCond_Intl;   
   BoundaryValues[0] = BoundValue;
   BoundaryValues[1] = BoundValue_Intl;

   TDiscreteForm3D *DiscreteFormMRhs_Galerkin;
   TDiscreteForm3D *DiscreteFormARhs_Galerkin; 
 //   TDiscreteForm3D *DiscreteFormMRhs_SUPG;
 //   TDiscreteForm3D *DiscreteFormARhs_SUPG;

  
   InitializeDiscreteFormsScalar(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormRhs, BilinearCoeffs);
  
     switch(Disctype)
      {
       case GALERKIN:
       case LOCAL_PROJECTION:
            DiscreteFormARhs = DiscreteFormARhs_Galerkin;
            DiscreteFormMRhs = DiscreteFormMRhs_Galerkin;
       break;
      
 //       case SUPG:
 //            DiscreteFormARhs = DiscreteFormARhs_SUPG;
 //            DiscreteFormMRhs = DiscreteFormMRhs_SUPG;
 //       break;
      
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
    
      RHSs[0] = RhsArray[i];
  
      // A Matrix     
      SQMATRICES[0] = sqmatrixA[i];
      AMatRhsAssemble[i] = new TAssembleMat3D(1, fesp, 1, SQMATRICES, 0, NULL, 1, RHSs, ferhs, 
                               DiscreteFormARhs, BoundaryConditions, BoundaryValues, aux);
      AMatRhsAssemble[i]->Init();
  
      // M matrix     
      SQMATRICES[0] = sqmatrixM[i];
      MMatRhsAssemble[i] = new TAssembleMat3D(1, fesp, 1, SQMATRICES, 0, NULL, 1, RHSs, ferhs, 
                               DiscreteFormMRhs, BoundaryConditions, BoundaryValues, aux);
      MMatRhsAssemble[i]->Init();   
      
      //setup the multigrid solver
      if(SOLVER==GMG)
       {
 #ifdef _MPI  
        MGLevel = new TMGLevel3D(i, SQMATRICES[0], RHSs[0], SolArray[i], ParComm[i], ParMapper[i], N_aux, NULL);
 #else
        MGLevel = new TMGLevel3D(i, SQMATRICES[0], RHSs[0], SolArray[i], N_aux, NULL);
 #endif
        MG->AddLevel(MGLevel);
       }  

     } // for(i=Star 

/** by default the intenal space is solved at the "N_Levels"th level of the physical space */
    PhycicalSpace_Level = N_Levels-1;

  // if(TDatabase::ParamDB->REACTOR_P5==0)
    // {
    // inernal eqn will be solved at every dof of space
    N_PhysicalPts = FeSpaces[PhycicalSpace_Level]->GetN_DegreesOfFreedom();
    PhysicalPts = new double[3*N_PhysicalPts];
    PhysicalPos = new TFEVectFunct3D(FeSpaces[PhycicalSpace_Level], PosString, PosString, PhysicalPts, N_PhysicalPts, 3);

    PhysicalPos->GridToData();
    Physical_X = PhysicalPts;
    Physical_Y = PhysicalPts+N_PhysicalPts;
    Physical_Z = PhysicalPts+N_PhysicalPts;           
    // }
  //  else
    // {
    //   GetPhysicalSpaceQuadPts(N_PhysicalPts, Scalar_Spaces[PBE_INDEX], IntlX, IntlY, IntlZ);  
    // }

   cout << " N_PhySpacePts : " << N_PhysicalPts <<endl;
  //=========================================================================
  /** Init for internal spaces */
  //=========================================================================  

    /** set data for multigrid */
  mg_type_intl = TDatabase::ParamDB->SC_MG_TYPE_SCALAR_INTL;
 
  if(TDatabase::ParamDB->SOLVER_TYPE_INTL==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE_INTL==DIRECT)
   { 
     mg_type_intl=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SCALAR_INTL = mg_type_intl;
    }
  
  if(mg_type_intl)
   {
    mg_level_intl =  LEVELS_INTL + 1;
    ORDER_INTL = -1;
   }
  else
   {
    mg_level_intl = LEVELS_INTL;
    ORDER_INTL = TDatabase::ParamDB->ANSATZ_ORDER_INTL;
   }

  FeSpaces_Intl = new TFESpace3D*[mg_level_intl];  
  FeFunctions_Intl = new TFEFunction3D*[mg_level_intl]; 
  Sol_array_Intl = new double*[mg_level_intl];
  Rhs_array_Intl = new double*[mg_level_intl];

//=========================================================================
// construct all internal finite element spaces
//=========================================================================
  for(i=0;i<LEVELS_INTL;i++)
   {   
    if(i)
     { Domain_Intl->RegRefineAll(); }
     
#ifdef _MPI
     if(rank == out_rank)
       printf("Level :: %d\n\n",i);
     if(i)
     {
       Domain_Intl->GenerateEdgeInfo();
       Domain_Crop(Comm, Domain_Intl);  // removing unwanted cells in the hallo after refinement 
     }
#endif

     Coll_Intl=Domain_Intl->GetCollection(It_Finest, 0);
     
     // fespaces for internal equation     
     FeSpaces_Intl[i] =  new TFESpace3D(Coll_Intl, Name, Description, BoundCond_Intl, ORDER_INTL);     

#ifdef _MPI
     FeSpaces_Intl[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
     
     //multilevel multigrid disc
     if(i==LEVELS_INTL-1 && mg_type_intl==1) 
      {
       ORDER_INTL = TDatabase::ParamDB->ANSATZ_ORDER_INTL;
       FeSpaces_Intl[mg_level_intl-1] =  new TFESpace3D(Coll_Intl, Name, Description, BoundCond_Intl, ORDER_INTL);
       if(ORDER_INTL==0)
         FeSpaces_Intl[mg_level_intl-1]->SetAsDGSpace();
#ifdef _MPI
       FeSpaces_Intl[mg_level_intl-1]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
      } //  if(i==LEVELS-1 && i!=mg_level_intl-1) 

if(ORDER_INTL==0)
 FeSpaces_Intl[i]->SetAsDGSpace();

    N_DoF_Intl = FeSpaces_Intl[i]->GetN_DegreesOfFreedom();
    sol_intl = new double[N_DoF_Intl];
    rhs_intl = new double[N_DoF_Intl];
    Sol_array_Intl[i] = sol_intl;
    Rhs_array_Intl[i] = rhs_intl;   


    FeFunction_Intl  = new TFEFunction3D(FeSpaces_Intl[i], Cstring, Cstring, sol_intl, N_DoF_Intl);  
    FeFunctions_Intl[i] = FeFunction_Intl;
     
    if(i==LEVELS_INTL-1 && mg_type_intl==1) 
     {  
      N_DoF_Intl = FeSpaces_Intl[mg_level_intl-1]->GetN_DegreesOfFreedom();
      sol_intl = new double[N_DoF_Intl];
      rhs_intl = new double[N_DoF_Intl];
      Sol_array_Intl[mg_level_intl-1] = sol_intl;
      Rhs_array_Intl[mg_level_intl-1] = rhs_intl;  

      FeFunction_Intl  = new TFEFunction3D(FeSpaces_Intl[i], Cstring, Cstring, sol_intl, N_DoF_Intl);   
      FeFunctions_Intl[mg_level_intl-1] = FeFunction_Intl;
     }//   if(i==LEVELS_INTL-1 && mg_type_intl==1) 

#ifdef _MPI
     N_Cells = Coll_Intl->GetN_Cells();
     printf("rank=%d\t N_Cells   : %d\t Dof all   :%d\n",rank,N_Cells,N_DOF);
#endif
   }// for(i=0;i<LEVELS_INTL;i++)

/** by default the physcal space is solved at the "LEVELS_INTL"th level of the internal space */
  InternalSpace_Level = LEVELS_INTL-1;

  N_InternalPts = FeSpaces_Intl[InternalSpace_Level]->GetN_DegreesOfFreedom();
  InternalPts = new double[3*N_InternalPts];
  InternalPos = new TFEVectFunct3D(FeSpaces_Intl[InternalSpace_Level], PosString, PosString, InternalPts, N_InternalPts, 3);

  InternalPos->GridToData();
  Internal_X = InternalPts;
  Internal_Y = InternalPts+N_InternalPts;
  Internal_Z = InternalPts+N_InternalPts;   

   cout << " N_InternalPts : " << N_InternalPts <<endl;


   cout << "SystemTCD3D_3C.C internal Init done! " <<endl;
   // beyod this point, yet to be implemented 
   exit(0);
} // Init


// void GetPhysicalSpaceQuadPts(int &N_PhySpacePts, TFESpace3D *FESpace3D, double *&IntlX, double *&IntlY, double *&IntlZ)
// {
//   int i,j,k,l, m;
//   int N_Cells, PolynomialDegree, N_Points;

//   double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D], AbsDetjk[MaxN_QuadPoints_3D];
//   double *xi, *eta, *zeta, *weights;

//   TBaseCell *cell;
//   TCollection *Coll;
//   FE3D FEId;
//   TFE3D *Element;
//   BF3DRefElements RefElement;
//   QuadFormula3D QuadFormula;
//   TQuadFormula3D *qf3;
//   TRefTrans3D *F_K;

//   Coll = FESpace3D->GetCollection();
//   N_Cells = Coll->GetN_Cells();

//    m = 0;
//   for(i=0; i<N_Cells; i++)
//   {
//     cell = Coll->GetCell(i);
//     FEId = FESpace3D->GetFE3D(i, cell);
//     Element = TFEDatabase3D::GetFE3D(FEId);
//     RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
//     PolynomialDegree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);

//      switch(RefElement)
//     {
//       case BFUnitHexahedron:
//         QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(3*PolynomialDegree);
//         F_K = TFEDatabase3D::GetRefTrans3D(HexaAffin);
//         ((THexaAffin *)F_K)->SetCell(cell);
//         break;

//       case BFUnitTetrahedron:
//         QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(3*PolynomialDegree-1);
//         F_K = TFEDatabase3D::GetRefTrans3D(TetraAffin);
//         ((TTetraAffin *)F_K)->SetCell(cell);
//         break;
//     }                                             // endswitch

// //     cout << "QuadFormula: " << QuadFormula << endl;
//     qf3 = TFEDatabase3D::GetQuadFormula3D(QuadFormula);
//     qf3->GetFormulaData(N_Points, weights, xi, eta, zeta);

//     if(i==0)
//     {
//       N_PhySpacePts = N_Points*N_Cells;

//       IntlX = new double[N_PhySpacePts];
//       IntlY = new double[N_PhySpacePts];
//       IntlZ = new double[N_PhySpacePts];      
//     }

//     switch(RefElement)
//     {
//       case BFUnitHexahedron:
//         ((THexaAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, AbsDetjk);
//         break;

//       case BFUnitTetrahedron:
//         ((TTetraAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, AbsDetjk);
//         break;
//     }     

//     for(j=0; j<N_Points; j++)
//     {
//       IntlX[m] = X[j];
//       IntlY[m] = Y[j];
//       IntlZ[m] = Z[j];
//       m++;
//     }
       
//   } //   for(i=0; i<N_Cells; i++)
  
// //      cout<< N_PhySpacePts << " N_PhySpacePts " << m <<endl;
// //     for(i=0; i<N_PhySpacePts; i++)
// //      cout<< i << " IntlX " << IntlX[i] << " IntlY " << IntlY[i]  << " IntlZ " << IntlZ[i]  <<  endl;
// //    exit(0);
// //   cout << " GetPhysicalSpaceQuadPts " << endl;
// //   exit(0);
// } // GetPhysicalSpaceQuadPts



// void TSystemTCD3D::AssembleMRhs()
// {
//   //this is set to true for direct solver factorization
//   factorize = true;
  
//   int i, N_DOF_low, N_Active;

//    for(i=Start_Level;i<N_Levels;i++)
//     {      
//      N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
//      N_Active =  FeSpaces[i]->GetActiveBound();
  
//      /** initialize matrices and rhs */
//      MMatRhsAssemble[i]->Reset(); 
     
//      // assemble
//      MMatRhsAssemble[i]->Assemble3D();
     
//      /** free the Mass mat array, no need in time loop */
//      MMatRhsAssemble[i]->DeAllocate();
     
//      /** set rhs for Dirichlet nodes */
//      memcpy(SolArray[i]+N_Active, RhsArray[i]+N_Active, (N_DOF_low - N_Active)*SizeOfDouble);         
//     } //  for(i=Start_Level;i<N_Levels;i++)

// } // TSystemMatScalar3D::AssembleMRhs 


// void TSystemTCD3D::AssembleARhs()
// {
//   //this is set to true for direct solver factorization
//   factorize = true;
  
//   int i, N_DOF_low, N_Active;

//    for(i=Start_Level;i<N_Levels;i++)
//     {    
//      N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
//      N_Active =  FeSpaces[i]->GetActiveBound();

//      /** reset the matrix and rhs */
//      AMatRhsAssemble[i]->Reset(); 
    
//      // assemble
//      AMatRhsAssemble[i]->Assemble3D();     

//      /** set rhs for Dirichlet nodes */
//      memcpy(SolArray[i]+N_Active, RhsArray[i]+N_Active, (N_DOF_low - N_Active)*SizeOfDouble);           
//     }//   for(i=Start_Level;i<N_Le  
    
// } // TSystemMatScalar3D::AssembleARhs 

// void TSystemTCD3D::AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol
// #ifdef _MPI
//                                              , double **Rhs_array
// #endif
//                                              )
// {
//     int i, N_Active;
//     double tau;
    
//     if(SystMatAssembled)
//      {
//       OutPut("System is has to be restored before calling AssembleSystMat! " <<endl);
//       exit(0);
//      }
    
//     SQMATRICES[0] = sqmatrixM[N_Levels-1];

//     N_Active =  FeSpaces[N_Levels-1]->GetActiveBound();     
//     tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;  
    
//     memset(B, 0, N_DOF*SizeOfDouble); 

//     /** old rhs multiplied with current subtime step and theta3 on B */
//     Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3,  oldrhs, B);    

//     /** add rhs from current sub time step to rhs array B */
//     Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4,  rhs, B);    

//     /** M = M + (- tau*THETA2)A */
//      MatAdd(sqmatrixM[N_Levels-1], sqmatrixA[N_Levels-1], - tau*TDatabase::TimeDB->THETA2);
//      gamma = -tau*TDatabase::TimeDB->THETA2;  // set current factor of steady state matrix

//      /** defect = M * oldsol */
//      memset(defect, 0, N_DOF*SizeOfDouble);  
//      MatVectActive(sqmatrixM[N_Levels-1], oldsol, defect); 
//     //cout << "defect " << Ddot(N_Active, sol, sol)<< endl;      
    
//      /** B:= B + defec  */
//      Daxpy(N_Active, 1, defect, B);

//      /** set Dirichlet values */
//      memcpy(B+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);  
//      memcpy(sol+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
          
//      /** assemble the system matrix */
//      for(i=Start_Level;i<N_Levels;i++)   
//       {
//        if(i==N_Levels-1)
//          { MatAdd(sqmatrixM[i], sqmatrixA[i], -gamma + tau*TDatabase::TimeDB->THETA1);}
//         else
//          { MatAdd(sqmatrixM[i], sqmatrixA[i], tau*TDatabase::TimeDB->THETA1);} 
         
// #ifdef _MPI  
//        SQMATRICES[0] = sqmatrixM[i];  
// #endif
//        }
//      gamma = tau*TDatabase::TimeDB->THETA1;
     
// //have to shift this in pardirectsolver     
// #ifdef _OMPONLY     
//     if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
//       DS->AssembleMatrix(sqmatrixM[N_Levels-1]);
// #endif
     
//      SystMatAssembled  = TRUE;

// } // AssembleSystMat

// void TSystemTCD3D::RestoreMassMat()
// {
//  int i;

//   if(SystMatAssembled)
//    {
//      // restore the mass matrix
//      for(i=Start_Level;i<N_Levels;i++)  
//       MatAdd(sqmatrixM[i], sqmatrixA[i], -gamma);
     
//      gamma = 0.;
//      SystMatAssembled  = FALSE;
//    }
//   else
//   {
//     cout << "System is not assembled to restore " <<endl;
//     exit(0);
//   }

// }

// void TSystemTCD3D::Solve(double *sol)
// {  
//     switch(SOLVER)
//      {
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

//       case DIRECT:
// #ifdef _MPI
// 	TDS->Solve(sol, B, factorize);
// // 	exit(0);
// #endif

// #ifdef _OMPONLY
// 	if(TDatabase::ParamDB->DSType == 1)
// 	  DS->Solve(sol, B, factorize);
// 	else{
// 	  OutPut("Select Proper Solver" << endl);
// 	  exit(0);
// 	}
// #endif

// #ifdef _SEQ
//         DirectSolver(sqmatrixM[N_Levels-1], B, sol);
// #endif
// 	//this is set to false for direct solver factorization
//         factorize = false;
//       break;      
 
//       default:
//             OutPut("Unknown Solver" << endl);
//             exit(4711);;
//      }
     
// }

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



