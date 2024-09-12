/** ************************************************************************
 * @brief     source file for TSystemTCD3D_3DSuraface
 * @author    Sashikumaar Ganesan
 * @date      07.10.2019
 * @History
 ************************************************************************  */
#include <Database.h>
#include <SystemPBE3D.h>
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

// Added for the 1D Mesh Generation 
#include <MacroCell.h>
#include <Joint.h>
#include <JointEqN.h>

// Added for internal System
#include <SystemADI1D.h>

TSystemPBE1D::TSystemPBE1D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver,
                           TFEFunction3D **spatialFEFunctions, DoubleFunctND *getKernel, DoubleFunctND *getRhs, TDomain *Domain_Intl)
    : TSystemTCD3D(N_levels, fespaces, sol, rhs, disctype, solver)
{
  char *GEO, *PRM;
  int N_Cells, N_dof_level;
#ifdef _SMPI
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  MPI_Comm_size(TDatabase::ParamDB->Comm, &mpi_size);
#endif

  // Thivin - Obtain the externally generated - internal Mesh co-ordinates
  Domain_Intl = new TDomain();
  TCollection *Coll_Intl;
  int N_Cells_Intl;
  double L0, L1;
  double *IntlPosL;

  // Spatial position for 1D Data for Internal Co-ordinates
  int N = 5;
  L0 = 0.1;
  L1 = 1.0;

  Generate1DMesh(Domain_Intl, L0, L1, N);

  Coll_Intl = Domain_Intl->GetCollection(It_Finest, 0);
  N_Cells_Intl = Coll_Intl->GetN_Cells();

  // create an ADISystem1D object for the internal domain


  // TSystemADI1D* ADISystem1D = new TSystemADI1D(N_levels, L0, L1, );     


  B1 = new TSquareMatrix3D *[N_Levels];
  B2 = new TSquareMatrix3D *[N_Levels];
  B3 = new TSquareMatrix3D *[N_Levels];

  if (Disctype == SUPG)
  {
    S1_SUPG = new TSquareMatrix3D *[N_Levels];
    S2_SUPG = new TSquareMatrix3D *[N_Levels];
    S3_SUPG = new TSquareMatrix3D *[N_Levels];
    S1S2_SUPG = new TSquareMatrix3D *[N_Levels];
    S1S3_SUPG = new TSquareMatrix3D *[N_Levels];
    S2S3_SUPG = new TSquareMatrix3D *[N_Levels];
    M1_SUPG = new TSquareMatrix3D *[N_Levels];
    M2_SUPG = new TSquareMatrix3D *[N_Levels];
    M3_SUPG = new TSquareMatrix3D *[N_Levels];
  }

  for (int i = Start_Level; i < N_Levels; i++)
  {
    B1[i] = new TSquareMatrix3D(sqstructure[i]);
    B2[i] = new TSquareMatrix3D(sqstructure[i]);
    B3[i] = new TSquareMatrix3D(sqstructure[i]);

    if (Disctype == SUPG)
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

  if (Disctype == SUPG)
    for (int i = 0; i < 9; i++)
      Rhs_SUPG[i] = new double *[N_Levels];

  for (int i = Start_Level; i < N_Levels; i++)
  {
    N_dof_level = FeSpaces[i]->GetN_DegreesOfFreedom();
    Rhs1Array[i] = new double[N_dof_level];
    Rhs2Array[i] = new double[N_dof_level];
    Rhs3Array[i] = new double[N_dof_level];
  }

  if (Disctype == SUPG)
  {
    for (int j = 0; j < 9; j++)
      for (int i = Start_Level; i < N_Levels; i++)
        Rhs_SUPG[j][i] = new double[N_dof_level];
  }

  rhs1old = new double[N_DOF];
  rhs2old = new double[N_DOF];
  rhs3old = new double[N_DOF];

  RhsAssemble = new TAssembleMat3D *[N_Levels];

  /** need it for solver */
  sqmatrices = (TSquareMatrix **)SQMATRICES_RTE;

  SpatialFEFunctions = spatialFEFunctions;

  GetKernel = getKernel;
  GetRhs = getRhs;
  // cout << "Domain_Intl " << N_SurfCells <<endl;

  // exit(0);

} // constructor

TSystemPBE1D::~TSystemPBE1D()
{
  int i;

  for (i = Start_Level; i < N_Levels; i++)
  {
    delete sqstructure[i];
    delete sqmatrixA[i];
  }

  delete[] sqstructure;
  delete[] sqmatrixA;

  if (SOLVER == GMG && TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
  {
    delete[] Itmethod_sol;
    delete[] Itmethod_rhs;
  }

  delete[] B;
  delete[] defect;

  //   delete sqmatrixM;
  //   if(Disctype==SDFEM || Disctype==SUPG)
  //    {
  //      delete sqmatrixS;
  //      delete sqmatrixK;
  //    }
}

void TSystemPBE1D::Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue,
                        BoundCondFunct3D *BoundCond_Intl, BoundValueFunct3D *BoundValue_Intl, TAuxParam3D *aux,
                        BoundCondFunct2D *SurfBDCond, BoundValueFunct2D *SurfBDValue)
{

  char Name[] = "IntenalSystem";
  char Description[] = "InternalSystem";
  char Cstring[] = "U";
  char PosString[] = "PhysicalSpaceCoord";

#ifdef _MPI
  if (SOLVER == DIRECT)
  {
    SQMATRICES_RTE[0] = sqmatrixM[N_Levels - 1];
    TDS = new TParDirectSolver(ParComm[N_Levels - 1], NULL, SQMATRICES_RTE, NULL);
  }
#endif

#ifdef _OMPONLY
  if (SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
  {
    DS = new TParDirectSolver(sqmatrixM[N_Levels - 1]);
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

  TDiscreteForm3D *DiscreteFormMRhs_Galerkin;
  TDiscreteForm3D *DiscreteFormARhs_Galerkin;
  //   TDiscreteForm3D *DiscreteFormMRhs_SUPG;
  //   TDiscreteForm3D *DiscreteFormARhs_SUPG;

  InitializeDiscreteFormsScalar(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormRhs, BilinearCoeffs);

  switch (Disctype)
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
    exit(4711);
    ;
  }

  // initialize the assemble
  if (aux == NULL)
  {
    aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
  }

  for (i = Start_Level; i < N_Levels; i++)
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
    RhsAssemble[i] = new TAssembleMat3D(1, fesp, 1, SQMATRICES, 0, NULL, 1, RHSs, ferhs,
                                        DiscreteFormMRhs, BoundaryConditions, BoundaryValues, aux);
    RhsAssemble[i]->Init();

    // setup the multigrid solver
    if (SOLVER == GMG)
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
  SpatialSpace_Level = N_Levels - 1;

  N_SpatialPts = FeSpaces[SpatialSpace_Level]->GetN_DegreesOfFreedom();
  SpatialPts = new double[3 * N_SpatialPts];
  SpatialPos = new TFEVectFunct3D(FeSpaces[SpatialSpace_Level], PosString, PosString, SpatialPts, N_SpatialPts, 3);

  SpatialPos->GridToData();
  Spatial_X = SpatialPts;
  Spatial_Y = SpatialPts + N_SpatialPts;
  Spatial_Z = SpatialPts + N_SpatialPts;

  //  cout << " N_SpatialPts : " << N_SpatialPts <<endl;
  //=========================================================================
  /** Init for internal spaces */
  //=========================================================================

  /** set data for multigrid */
  ORDER_INTL = TDatabase::ParamDB->ANSATZ_ORDER_INTL;

  //===================================================== ====================
  // construct all internal finite element spaces
  //=========================================================================
  // fe (bulk) space for internal equation, needed for interpolation and output
  FeSpaces_Intl1D = new TFESpace1D(Coll_Intl, "Internal-1D Space", "Internal-1D Space", TDatabase::ParamDB->TEST_ORDER);
  N_DoF_Intl1D = FeSpaces_Intl1D->GetN_DegreesOfFreedom();
  Sol_Intl1D = new double[N_DoF_Intl1D];

  FeFunction_Intl1D = new TFEFunction1D(FeSpaces_Intl1D, Cstring, Cstring, Sol_Intl1D, N_DoF_Intl1D);

  // std::ostringstream os;
  // os << " ";
  // Output = new TOutput1D(1, 1, 1, 1, Domain_Intl);
  // os.seekp(std::ios::beg);
  // Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
  // Output->AddFEFunction(FeFunction_Intl1D);

  // os.seekp(std::ios::beg);
  // os << TDatabase::ParamDB->VTKBASENAME<<".vtk" << ends;
  // Output->WriteVtk(os.str().c_str());

  // ---- THIVIN -- Commented out below ------------------------------ //

  // fe (surface) space for internal equation
  /** physcal space is solved on the surface of the internal domain */
  // cout << "Number of Surface Cells: " << N_SurfCells<<endl;

  // FeSpaces_Intl = new TFESpace2D(Surf_Coll, Name, Description, SurfBDCond, ORDER_INTL, NULL);
  // N_DoF_Intl = FeSpaces_Intl->GetN_DegreesOfFreedom();
  // cout << "N_Surf_Dof : " << N_DoF_Intl <<endl;
  // cout << "ORDER_INTL : " << ORDER_INTL <<endl;

  // sol_intl = new double[N_DoF_Intl*N_SpatialPts];
  // rhs_intl = new double[N_DoF_Intl*N_SpatialPts];
  // BE_MatScale = new double[N_DoF_Intl];
  // //get the coordinates of the internal points (surface points), needed for physical space convection term
  // InternalPts = new double[3*N_DoF_Intl];

  // FeSpaces_Intl1D->GetSurfDOFPosition(N_SurfCells, Cell_IDX, Joint_IDX, InternalPts);

  //  cout << " N_DoF_Intl : " << N_DoF_Intl <<endl;
  //  for(i=0;i<N_DoF_Intl;i++)
  //  cout << " InternalPts : " << InternalPts[3*i] << " " << InternalPts[3*i+1]<< " " << InternalPts[3*i+2]<<endl;

  // cout << "SystemTCD3D_3C.C internal Init done! " <<endl;

} // Init

void TSystemPBE1D::Interpolate(DoubleFunct3D *Exact, double *Sol)
{
  double *sol = SpatialFEFunctions[SpatialSpace_Level]->GetValues();

  for (int i = 0; i < N_DoF_Intl; ++i)
  {
    TDatabase::ParamDB->P11 = InternalPts[3 * i];
    TDatabase::ParamDB->P12 = InternalPts[3 * i + 1];
    TDatabase::ParamDB->P13 = InternalPts[3 * i + 2];
    SpatialFEFunctions[SpatialSpace_Level]->Interpolate(Exact);

    memcpy(Sol + (i * N_SpatialPts), sol, N_SpatialPts * SizeOfDouble);
  }
}

void TSystemPBE1D::AssembleABRhs()
{
  // this is set to true for direct solver factorization
  factorize = true;

  int i, N_DOF_low, N_Active;

  for (i = Start_Level; i < N_Levels; i++)
  {
    /** reset the matrix and rhs */
    AMatRhsAssemble[i]->Reset();

    // assemble
    AMatRhsAssemble[i]->Assemble3D();

  } //   for(i=Start_Level;i<N_Le

} // TSystemMatScalar3D::AssembleARhs

void TSystemPBE1D::AssembleRhs()
{
  int i, N_DOF_low, N_Active;

  memcpy(rhs1old, Rhs1Array[N_Levels - 1], N_DOF * SizeOfDouble);
  memcpy(rhs2old, Rhs2Array[N_Levels - 1], N_DOF * SizeOfDouble);
  memcpy(rhs3old, Rhs3Array[N_Levels - 1], N_DOF * SizeOfDouble);

  for (i = Start_Level; i < N_Levels; i++)
  {
    /** reset the matrix and rhs */
    RhsAssemble[i]->Reset();

    // assemble
    RhsAssemble[i]->Assemble3D();
  } //   for(i=Start_Level;i<N_Le

} // AssembleRhs()

void TSystemPBE1D::AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol,
                                   double b1, double b2, double b3, CoeffFct3D *SystemRhs
#ifdef _MPI
                                   ,
                                   double **Rhs_array
#endif
)
{
  int i, N_Active;
  double tau, param[3];
  // double diffusion = TDatabase::ParamDB->PE_NR;

  if (SystMatAssembled)
  {
    OutPut("System has to be restored before calling AssembleSystMat! " << endl);
    exit(0);
  }

  N_Active = FeSpaces[N_Levels - 1]->GetActiveBound();
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  memset(B, 0, N_DOF * SizeOfDouble);

  /** old rhs multiplied with current subtime step and theta3 on B */
  if (fabs(TDatabase::TimeDB->THETA3) > 0)
  {
    Daxpy(N_Active, tau * TDatabase::TimeDB->THETA3, oldrhs, B);
    Daxpy(N_Active, b1 * tau * TDatabase::TimeDB->THETA3, rhs1old, B);
    Daxpy(N_Active, b2 * tau * TDatabase::TimeDB->THETA3, rhs2old, B);
    Daxpy(N_Active, b3 * tau * TDatabase::TimeDB->THETA3, rhs3old, B);

    cout << "Old Rhs..Array in AssembleRhs() need to be stored befrore using this time-disc! " << endl;
    exit(0);
  }

  /** add rhs from current sub time step to rhs array B */
  if (fabs(TDatabase::TimeDB->THETA4) > 0)
  {
    RHSs[0] = rhs;
    RHSs[1] = Rhs1Array[N_Levels - 1];
    RHSs[2] = Rhs2Array[N_Levels - 1];
    RHSs[3] = Rhs3Array[N_Levels - 1];

    if (Disctype == SUPG)
    {
      for (int j = 0; j < 9; j++)
        RHSs[4 + j] = Rhs_SUPG[j][N_Levels - 1];
    }

    param[0] = b1;
    param[1] = b2;
    param[2] = b3;
    SystemRhs(N_Active, param, B, NULL, RHSs, NULL);
  }

  /** M + \delta_K (U, S\cdot\grad V) */
  if (Disctype == SUPG)
  {
    // time consistent SUPG
    MatAdd(sqmatrixM[N_Levels - 1], M1_SUPG[N_Levels - 1], b1);
    MatAdd(sqmatrixM[N_Levels - 1], M2_SUPG[N_Levels - 1], b2);
    MatAdd(sqmatrixM[N_Levels - 1], M3_SUPG[N_Levels - 1], b3);
  }

  /** M = M + (- tau*THETA2)A */
  if (fabs(TDatabase::TimeDB->THETA2) > 0)
  {
    MatAdd(sqmatrixM[N_Levels - 1], sqmatrixA[N_Levels - 1], -tau * TDatabase::TimeDB->THETA2);
    MatAdd(sqmatrixM[N_Levels - 1], B1[N_Levels - 1], -b1 * tau * TDatabase::TimeDB->THETA2);
    MatAdd(sqmatrixM[N_Levels - 1], B2[N_Levels - 1], -b2 * tau * TDatabase::TimeDB->THETA2);
    MatAdd(sqmatrixM[N_Levels - 1], B3[N_Levels - 1], -b3 * tau * TDatabase::TimeDB->THETA2);

    gammab0 = -tau * TDatabase::TimeDB->THETA2;      // set current factor of steady statmatrix
    gammab1 = -b1 * tau * TDatabase::TimeDB->THETA2; // set current factor of steady statmatrix
    gammab2 = -b2 * tau * TDatabase::TimeDB->THETA2; // set current factor of steady state matrix
    gammab3 = -b3 * tau * TDatabase::TimeDB->THETA2; // set current factor of steady state matrix

    if (Disctype == SUPG)
    {
      MatAdd(sqmatrixM[N_Levels - 1], S1_SUPG[N_Levels - 1], -b1 * b1 * tau * TDatabase::TimeDB->THETA2);
      MatAdd(sqmatrixM[N_Levels - 1], S2_SUPG[N_Levels - 1], -b2 * b2 * tau * TDatabase::TimeDB->THETA2);
      MatAdd(sqmatrixM[N_Levels - 1], S3_SUPG[N_Levels - 1], -b3 * b3 * tau * TDatabase::TimeDB->THETA2);

      MatAdd(sqmatrixM[N_Levels - 1], S1S2_SUPG[N_Levels - 1], -b1 * b2 * tau * TDatabase::TimeDB->THETA2);
      MatAdd(sqmatrixM[N_Levels - 1], S1S3_SUPG[N_Levels - 1], -b1 * b3 * tau * TDatabase::TimeDB->THETA2);
      MatAdd(sqmatrixM[N_Levels - 1], S2S3_SUPG[N_Levels - 1], -b2 * b3 * tau * TDatabase::TimeDB->THETA2);

      gammab4 = -b1 * b1 * tau * TDatabase::TimeDB->THETA2; // set current factor of steady statmatrix
      gammab5 = -b2 * b2 * tau * TDatabase::TimeDB->THETA2; // set current factor of steady state matrix
      gammab6 = -b3 * b3 * tau * TDatabase::TimeDB->THETA2; // set current factor of steady state matrix
      gammab7 = -b1 * b2 * tau * TDatabase::TimeDB->THETA2; // set current factor of steady statmatrix
      gammab8 = -b1 * b3 * tau * TDatabase::TimeDB->THETA2; // set current factor of steady state matrix
      gammab9 = -b2 * b3 * tau * TDatabase::TimeDB->THETA2; // set current factor of steady state matrix
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
  memset(defect, 0, N_DOF * SizeOfDouble);
  MatVectActive(sqmatrixM[N_Levels - 1], oldsol, defect);
  // cout << "defect " << Ddot(N_Active, sol, sol)<< endl;

  /** B:= B + defec  */
  Daxpy(N_Active, 1, defect, B);

  /** set Dirichlet values */
  memcpy(B + N_Active, rhs + N_Active, (N_DOF - N_Active) * SizeOfDouble);
  memcpy(sol + N_Active, rhs + N_Active, (N_DOF - N_Active) * SizeOfDouble);

  // for rteSpatialINternal.h example only (hard coded)
  Dscal((N_DOF - N_Active), b1 * b2 * b3, B + N_Active);
  Dscal((N_DOF - N_Active), b1 * b2 * b3, sol + N_Active);

  /** assemble the system matrix */
  for (i = Start_Level; i < N_Levels; i++)
  {
    if (i == N_Levels - 1)
    {

      MatAdd(sqmatrixM[i], sqmatrixA[i], -gammab0 + tau * TDatabase::TimeDB->THETA1);

      MatAdd(sqmatrixM[i], B1[i], -gammab1 + b1 * tau * TDatabase::TimeDB->THETA1);
      MatAdd(sqmatrixM[i], B2[i], -gammab2 + b2 * tau * TDatabase::TimeDB->THETA1);
      MatAdd(sqmatrixM[i], B3[i], -gammab3 + b3 * tau * TDatabase::TimeDB->THETA1);

      if (Disctype == SUPG)
      {
        MatAdd(sqmatrixM[i], S1_SUPG[i], -gammab4 + b1 * b1 * tau * TDatabase::TimeDB->THETA1);
        MatAdd(sqmatrixM[i], S2_SUPG[i], -gammab5 + b2 * b2 * tau * TDatabase::TimeDB->THETA1);
        MatAdd(sqmatrixM[i], S3_SUPG[i], -gammab6 + b3 * b3 * tau * TDatabase::TimeDB->THETA1);

        MatAdd(sqmatrixM[i], S1S2_SUPG[i], -gammab7 + b1 * b2 * tau * TDatabase::TimeDB->THETA1);
        MatAdd(sqmatrixM[i], S1S3_SUPG[i], -gammab8 + b1 * b3 * tau * TDatabase::TimeDB->THETA1);
        MatAdd(sqmatrixM[i], S2S3_SUPG[i], -gammab9 + b2 * b3 * tau * TDatabase::TimeDB->THETA1);
      }
    }
    else
    {
      MatAdd(sqmatrixM[i], sqmatrixA[i], tau * TDatabase::TimeDB->THETA1);

      MatAdd(sqmatrixM[i], B1[i], b1 * tau * TDatabase::TimeDB->THETA1);
      MatAdd(sqmatrixM[i], B2[i], b2 * tau * TDatabase::TimeDB->THETA1);
      MatAdd(sqmatrixM[i], B3[i], b3 * tau * TDatabase::TimeDB->THETA1);

      if (Disctype == SUPG)
      {
        // time consistent SUPG on low MG levels
        MatAdd(sqmatrixM[i], M1_SUPG[i], b1);
        MatAdd(sqmatrixM[i], M2_SUPG[i], b2);
        MatAdd(sqmatrixM[i], M3_SUPG[i], b3);

        MatAdd(sqmatrixM[i], S1_SUPG[i], b1 * b1 * tau * TDatabase::TimeDB->THETA1);
        MatAdd(sqmatrixM[i], S2_SUPG[i], b2 * b2 * tau * TDatabase::TimeDB->THETA1);
        MatAdd(sqmatrixM[i], S3_SUPG[i], b3 * b3 * tau * TDatabase::TimeDB->THETA1);

        MatAdd(sqmatrixM[i], S1S2_SUPG[i], b1 * b2 * tau * TDatabase::TimeDB->THETA1);
        MatAdd(sqmatrixM[i], S1S3_SUPG[i], b1 * b3 * tau * TDatabase::TimeDB->THETA1);
        MatAdd(sqmatrixM[i], S2S3_SUPG[i], b2 * b3 * tau * TDatabase::TimeDB->THETA1);
      }
    }

#ifdef _MPI
    SQMATRICES_RTE[0] = sqmatrixM[i];
#endif
  } // for(i=Start_Leve

  gammab0 = tau * TDatabase::TimeDB->THETA1;
  gammab1 = b1 * tau * TDatabase::TimeDB->THETA1; // set current factor of steady state matrix
  gammab2 = b2 * tau * TDatabase::TimeDB->THETA1; // set current factor of steady state matrix
  gammab3 = b3 * tau * TDatabase::TimeDB->THETA1; // set current factor of steady state matrix

  if (Disctype == SUPG)
  {
    gammab4 = b1 * b1 * tau * TDatabase::TimeDB->THETA1;
    gammab5 = b2 * b2 * tau * TDatabase::TimeDB->THETA1;
    gammab6 = b3 * b3 * tau * TDatabase::TimeDB->THETA1;
    gammab7 = b1 * b2 * tau * TDatabase::TimeDB->THETA1;
    gammab8 = b1 * b3 * tau * TDatabase::TimeDB->THETA1;
    gammab9 = b2 * b3 * tau * TDatabase::TimeDB->THETA1;
    gammab10 = b1;
    gammab11 = b2;
    gammab12 = b3;
  }

  // //have to shift this in pardirectsolver
  // #ifdef _OMPONLY
  //     if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
  //       DS->AssembleMatrix(sqmatrixM[N_Levels-1]);
  // #endif

  SystMatAssembled = TRUE;
  // cout << "System is assembled "  << b1 <<endl;
} // AssembleSystMat

void TSystemPBE1D::RestoreMassMat()
{
  int i;

  if (SystMatAssembled)
  {
    // restore the mass matrix
    for (i = Start_Level; i < N_Levels; i++)
    {
      MatAdd(sqmatrixM[i], sqmatrixA[i], -gammab0);

      MatAdd(sqmatrixM[i], B1[i], -gammab1);
      MatAdd(sqmatrixM[i], B2[i], -gammab2);
      MatAdd(sqmatrixM[i], B3[i], -gammab3);

      if (Disctype == SUPG)
      {

        // time consistent SUPG
        MatAdd(sqmatrixM[i], M1_SUPG[i], -gammab10);
        MatAdd(sqmatrixM[i], M2_SUPG[i], -gammab11);
        MatAdd(sqmatrixM[i], M3_SUPG[i], -gammab12);

        MatAdd(sqmatrixM[i], S1_SUPG[i], -gammab4);
        MatAdd(sqmatrixM[i], S2_SUPG[i], -gammab5);
        MatAdd(sqmatrixM[i], S3_SUPG[i], -gammab6);

        MatAdd(sqmatrixM[i], S1S2_SUPG[i], -gammab7);
        MatAdd(sqmatrixM[i], S1S3_SUPG[i], -gammab8);
        MatAdd(sqmatrixM[i], S2S3_SUPG[i], -gammab9);
      }
    }

    gammab0 = 0;
    gammab1 = 0.;
    gammab2 = 0.;
    gammab3 = 0.;
    if (Disctype == SUPG)
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

    SystMatAssembled = FALSE;
  }
  else
  {
    cout << "System is not assembled to restore " << endl;
    exit(0);
  }
}

void TSystemPBE1D::Solve(double *sol)
{
  switch (SOLVER)
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
    if (TDatabase::ParamDB->DSType == 1)
      DS->Solve(sol, B, factorize);
    else
    {
      OutPut("Select Proper Solver" << endl);
      exit(0);
    }
#endif

#ifdef _SEQ
    DirectSolver(sqmatrixM[N_Levels - 1], B, sol);

    // cout << "seq direct solver " <<endl;
#endif
    // this is set to false for direct solver factorization
    factorize = false;
    break;

  default:
    OutPut("Unknown Solver" << endl);
    exit(4711);
    ;
  }

} // solve

void TSystemPBE1D::GetQuadIntegral(int i, int &N_Points, double *&weights, double *x, double *y, double *z,
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
  {
    Isoparametric = TRUE;
  }

  RefTrans = TetraAffin;
  F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
  ((TTetraAffin *)F_K)->SetCell(Cell);
  QuadFormula = Gauss3Tria;
  qf2d = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
  qf2d->GetFormulaData(N_Points, weights, p1, p2);

  // compute normals for all quad points
  for (k = 0; k < N_Points; ++k)
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
      eta = 1 - p1[k] - p2[k];
      zeta = p2[k];
      break;
    case 3:
      xi = 0;
      eta = p1[k];
      zeta = p2[k];
      break;
    }

    ((TTetraAffin *)F_K)->GetTangentVectors(Joint_No, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
    n1[k] = (a2 * b3 - a3 * b2);
    n2[k] = (a3 * b1 - a1 * b3);
    n3[k] = (a1 * b2 - a2 * b1);

    ((TTetraAffin *)F_K)->GetOrigFromRef(xi, eta, zeta, x[k], y[k], z[k]);
  } //  for (k=0;k<N_Points;++

} // void TSystemPBE1D::GetQuadIntegral(i

void TSystemPBE1D::AssembleIntlMat(int CellNo, double tau, double *XPos, double &fact, double *InternalScaling)
{
  int i, k, kk, l, m, Cell_No, Joint_No, N_Points;
  double *weights, *weightsk, p1, *p2, xi, eta, zeta;
  double a1, a2, a3, b1, b2, b3, sigma_a, sigma_s;
  double val, len, len1, surf_area = 0., surf_areak, stiff;
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
  memset(InternalScaling, 0, N_SpatialPts * SizeOfDouble);

  if (N_DoF_Intl != N_SurfCells)
  {
    cout << "Only Po is allowed for surface mesh" << endl;
    exit(0);
  }

  i = CellNo;
  //  for(i=0;i<N_SurfCells; i++)
  {
    this->GetQuadIntegral(i, N_Points, weights, x, y, z, n1, n2, n3);

    mass = 0.;
    stiff = 0.;
    for (k = 0; k < N_Points; ++k)
    {
      len = sqrt(n1[k] * n1[k] + n2[k] * n2[k] + n3[k] * n3[k]);
      len *= weights[k]; /** X * n1*/
      ;

      mass += len;
      stiff += sigma_a * len;

      Coords[0] = x[k];
      Coords[1] = y[k];
      Coords[2] = z[k];
      // compute rhs
      for (m = 0; m < N_SpatialPts; ++m)
      {
        Coords[6] = XPos[m]; // spatial x
        Coords[7] = XPos[N_SpatialPts + m];
        Coords[8] = XPos[2 * N_SpatialPts + m];

        //  get rhs
        GetRhs(6, Coords, InnVal);

        InternalScaling[m] += len * InnVal[0]; //  rhs +=len*InnVal[0];
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

      stiff -= sigma_s * len * val;
      // surf_area += len;
    } // for (int k=

    // rhs scaling factor in backwar Euler
    // BE_MatScale[i] = (1. - tau*TDatabase::TimeDB->THETA2*stiff)/ (1. + tau*TDatabase::TimeDB->THETA1*stiff);

    fact = 1. / (1. + tau * TDatabase::TimeDB->THETA1 * stiff);
    // only backward Euler
    for (m = 0; m < N_SpatialPts; ++m)
    {
      InternalScaling[m] *= tau * fact;
      //  OutPut(fact << " Mat: " <<  InternalScaling[m] <<endl);
    }

  } // for(i=0;i
  // tau*TDatabase::TimeDB->THETA1
  //  OutPut(4.*Pi<< " Surf Area: " << surf_area <<endl);
  //  exit(0);

} // void TSystemPBE1D::AssembleIntlM

void TSystemPBE1D::GetErrors(DoubleFunctVect *Exact, TFEFunction3D *ScalarFunction, double *sol_all, double &l2)
{

  // delete [] Sol_AllL;
  // delete [] InternalError;

} // void TSystemPBE1D::GetErrors

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

// Objective : Generate 1D uniform internal co-ordinates
void TSystemPBE1D::Generate1DMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, z;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell **CellTree;

  N_V = N_Cells + 1;

  Lines = new int[2 * N_Cells];
  Vetrex = new TVertex *[N_V];


  // Generate uniform values of X from Start to End for N_Cells
  double *X = new double[N_V];
  h = (End - Start) / N_Cells;
  for (i = 0; i < N_V; i++)
    X[i] = Start + i * h;

  // print the values of X
  // for (i = 0; i < N_V; i++)
  //   cout << " X " << X[i] << endl;
  // exit(0);


  y = 0.;
  z = 0.;
  for (i = 0; i < N_Cells; i++)
  {
    Lines[2 * i] = i;
    Lines[2 * i + 1] = i + 1;
    Vetrex[i] = new TVertex(X[i], y, z);
  }

  Vetrex[N_Cells] = new TVertex(X[N_V - 1], y, z);

  CellTree = new TBaseCell *[N_Cells];

  for (i = 0; i < N_Cells; i++)
  {
    //     Vetrex[ i ]->GetCoords(x, y);
    //     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[Lines[2 * i]]);
    CellTree[i]->SetVertex(1, Vetrex[Lines[2 * i + 1]]);
    ((TMacroCell *)CellTree[i])->SetSubGridID(0);

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

  for (i = 1; i < N_Cells; i++)
  {
    Joint = new TJointEqN(CellTree[i - 1], CellTree[i]);

    CellTree[i - 1]->SetJoint(1, Joint);
    CellTree[i]->SetJoint(0, Joint);
  } // for(i=0;i<N_Cells;i++)

  // end joint(vertex)
  Joint = new TJointEqN(CellTree[N_Cells - 1]);
  CellTree[N_Cells - 1]->SetJoint(1, Joint);

  delete[] Lines;
}
