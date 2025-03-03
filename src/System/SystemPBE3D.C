/** ************************************************************************
 * @brief     source file for TSystemPBE3D
 * @author    Thivin Anandh, Sashikumaar Ganesan
 * @date      24.01.15
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

#include <FESpace3D.h>
#include <FEVectFunct3D.h>

#include <stdlib.h>
#include <string.h>

#include <omp.h>

TSystemPBE3D::TSystemPBE3D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver)
    : TSystemCD3D(N_levels, fespaces, sol, rhs, disctype, solver)
{
  int i;

  /** M mass and system matrix */
  sqmatrixM = new TSquareMatrix3D *[N_Levels];
  N_Matrices++;

  MMatRhsAssemble = new TAssembleMat3D *[N_Levels];

  for (i = Start_Level; i < N_Levels; i++)
  {
    sqmatrixM[i] = new TSquareMatrix3D(sqstructure[i]);
  }

  /** working rhs, used in AssembleSystMat() */
  B = new double[N_DOF];
  defect = new double[N_DOF];

  gamma = 0.;

  //   /** time-consistent part of the SUPG matrix */
  //   if(Disctype==SDFEM || Disctype==SUPG)
  //    {
  //     sqmatrixS = new TSquareMatrix3D(sqstructure);
  //     N_Matrices++;
  //
  //     sqmatrixK = new TSquareMatrix3D(sqstructure);
  //     N_Matrices++;
  //    }

  SystMatAssembled = FALSE;

} // constructor

TSystemPBE3D::~TSystemPBE3D()
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

void TSystemPBE3D::Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue,
                        TAuxParam3D *aux)
{
#ifdef _MPI
  if (SOLVER == DIRECT)
  {
    SQMATRICES[0] = sqmatrixM[N_Levels - 1];
    TDS = new TParDirectSolver(ParComm[N_Levels - 1], NULL, SQMATRICES, NULL);
  }
#endif

#ifdef _OMPONLY
  if (SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
  {
    DS = new TParDirectSolver(sqmatrixM[N_Levels - 1]);
  }
#endif

  int i;

  BoundaryConditions[0] = BoundCond;
  BoundaryValues[0] = BoundValue;

  TDiscreteForm3D *DiscreteFormMRhs_Galerkin;
  TDiscreteForm3D *DiscreteFormARhs_Galerkin;
  //   TDiscreteForm3D *DiscreteFormMRhs_SUPG;
  //   TDiscreteForm3D *DiscreteFormARhs_SUPG;

  InitializeDiscreteFormsScalar(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormRhs, BilinearCoeffs);

  switch (Disctype)
  {
  case GALERKIN:
    //       case LOCAL_PROJECTION:
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
    MMatRhsAssemble[i] = new TAssembleMat3D(1, fesp, 1, SQMATRICES, 0, NULL, 1, RHSs, ferhs,
                                            DiscreteFormMRhs, BoundaryConditions, BoundaryValues, aux);
    MMatRhsAssemble[i]->Init();

    // setup the multigrid solver
    if (SOLVER == GMG)
    {
#ifdef _MPI
      MGLevel = new TMGLevel3D(i, SQMATRICES[0], RHSs[0], SolArray[i], ParComm[i], ParMapper[i], N_aux, NULL);
#else
      MGLevel = new TMGLevel3D(i, SQMATRICES[0], RHSs[0], SolArray[i], N_aux, NULL);
#endif
      MG->AddLevel(MGLevel);
    }

  } // for(i=Star
} // Init

// For Population Balance Equation with NSE Values.
// IT alsi includes the C . \grad(u_p)
void TSystemPBE3D::Init_with_NSEValues(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue,
                                       TAuxParam3D *aux)
{
  cout << "Correct Function for Init_with_NSEValues" << endl;
#ifdef _MPI
  if (SOLVER == DIRECT)
  {
    SQMATRICES[0] = sqmatrixM[N_Levels - 1];
    TDS = new TParDirectSolver(ParComm[N_Levels - 1], NULL, SQMATRICES, NULL);
  }
#endif

#ifdef _OMPONLY
  if (SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
  {
    DS = new TParDirectSolver(sqmatrixM[N_Levels - 1]);
  }
#endif

  int i;

  BoundaryConditions[0] = BoundCond;
  BoundaryValues[0] = BoundValue;

  TDiscreteForm3D *DiscreteFormMRhs_Galerkin;
  TDiscreteForm3D *DiscreteFormARhs_Galerkin;
  //   TDiscreteForm3D *DiscreteFormMRhs_SUPG;
  //   TDiscreteForm3D *DiscreteFormARhs_SUPG;

  InitializeDiscreteFormsScalarNSE(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormRhs, BilinearCoeffs);

  switch (Disctype)
  {
  case GALERKIN:
    //       case LOCAL_PROJECTION:
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
    MMatRhsAssemble[i] = new TAssembleMat3D(1, fesp, 1, SQMATRICES, 0, NULL, 1, RHSs, ferhs,
                                            DiscreteFormMRhs, BoundaryConditions, BoundaryValues, aux);
    MMatRhsAssemble[i]->Init();

    // setup the multigrid solver
    if (SOLVER == GMG)
    {
#ifdef _MPI
      MGLevel = new TMGLevel3D(i, SQMATRICES[0], RHSs[0], SolArray[i], ParComm[i], ParMapper[i], N_aux, NULL);
#else
      MGLevel = new TMGLevel3D(i, SQMATRICES[0], RHSs[0], SolArray[i], N_aux, NULL);
#endif
      MG->AddLevel(MGLevel);
    }

  } // for(i=Star
} // Init

// Function to import values for the drift velocity
void TSystemPBE3D::Init_for_drift_velocity(TFEVectFunct3D **fevect_drift_array, double *g, double fluid_rho, double *particle_rho, double fluid_viscosity,
                                           int N_internal, double *diameter_values, TFEVectFunct3D *fevect_b, int n_velocity_points)
{
  // Assign all the values
  m_fevect_drift_array = fevect_drift_array;
  m_g_array = g;
  m_fluid_rho = fluid_rho;
  m_particle_rho = particle_rho;
  m_fluid_viscosity = fluid_viscosity;
  m_N_internal = N_internal;
  m_diameter_values = diameter_values;
  m_fevect_fluid_velocity = fevect_b; // Stores the fevect function of the Underlying Fluid field.
  m_n_velocity_points = n_velocity_points;

  // Create a FeVectFunction to store the Intermediate Drift Velocity (v)
  // The actual drift velocity, which comes in stores the values of u_l, we are creating this to store the temporary values of v
  m_fevect_u_l = new TFEVectFunct3D *[m_N_internal];

  // Create an array to store the particle velocity similar to the drift velocity
  double **u_l = new double *[m_N_internal]();

  // Create with same FESpace as the drift velocity
  for (int i = 0; i < m_N_internal; i++)
  {
    u_l[i] = new double[3 * m_n_velocity_points]();

    m_fevect_u_l[i] = new TFEVectFunct3D(fevect_drift_array[i]->GetFESpace3D(), "u_particle", "u_particle", u_l[i], m_n_velocity_points, 3);
  }

  // Set up Array and memory to store the gradients of the drift velocity
  double **particle_velocity_gradient_x = new double *[m_N_internal];
  double **particle_velocity_gradient_y = new double *[m_N_internal];
  double **particle_velocity_gradient_z = new double *[m_N_internal];

  for (int i = 0; i < m_N_internal; i++)
  {
    particle_velocity_gradient_x[i] = new double[m_n_velocity_points]();
    particle_velocity_gradient_y[i] = new double[m_n_velocity_points]();
    particle_velocity_gradient_z[i] = new double[m_n_velocity_points]();
  }

  // obtain the co-ordinates of the Physical points where fluid velocity is stored.
  // For this get the FeVect of drift velocity for the first component and use GridToData to get the values
  // This is done to get the physical co-ordinates of the points where the fluid velocity is stored.
  m_fevect_drift = fevect_drift_array[0];
  m_physical_coordinates = new double[4 * m_n_velocity_points](); // here the 4th dimension will be used to store cell id.
  TFESpace3D *fespace_drift = m_fevect_drift->GetFESpace3D();

  // Create a new FeVectFunction to Store the Physical Co-ordinates
  TFEVectFunct3D *fevect_physical_coordinates = new TFEVectFunct3D(fespace_drift, "Physical_Coordinates", "Physical_Coordinates", m_physical_coordinates, m_n_velocity_points, 4);
  fevect_physical_coordinates->GridToDataWithCellid(); // populates cell id on 4th dimension
}

void TSystemPBE3D::SolveDriftVelocity(double timestep, int internal_level, double *solution)
{
  // Assume the incoming velocity is u_l ( which goes into PBE assembly)
  // We need to Solve for v(drift): which is u (fluid) + u_l (particle)
  // After solving for v(drift), we need to update the particle velocity u_l = v(drift) - u(fluid)
  // The equation to solve is d(v)/dt = -(v.grad)v - (1/\tau)(v - u) + (1 - mu)g

  // Obtain the drift velocity for the internal level
  TFEVectFunct3D *drift_velocity_fevect = m_fevect_drift_array[internal_level];

  // Store the Components of the drift velocity and the gradients
  TFEFunction3D *comp0 = drift_velocity_fevect->GetComponent(0);
  TFEFunction3D *comp1 = drift_velocity_fevect->GetComponent(1);
  TFEFunction3D *comp2 = drift_velocity_fevect->GetComponent(2);

  // Store the components
  TFEFunction3D *fluid_comp0 = m_fevect_fluid_velocity->GetComponent(0);
  TFEFunction3D *fluid_comp1 = m_fevect_fluid_velocity->GetComponent(1);
  TFEFunction3D *fluid_comp2 = m_fevect_fluid_velocity->GetComponent(2);

  // For ith Internal point, loop over all the physical points
  for (int index = 0; index < m_n_velocity_points; index++)
  {
    // obtain the gradients of the drift velocity at current point
    double x_coord = m_physical_coordinates[index];
    double y_coord = m_physical_coordinates[m_n_velocity_points + index];
    double z_coord = m_physical_coordinates[2 * m_n_velocity_points + index];
    int cell_id = (int)m_physical_coordinates[3 * m_n_velocity_points + index];
    TBaseCell *cell = m_fevect_drift->GetFESpace3D()->GetCollection()->GetCell(cell_id);

    cout << "-------------------------------------------- " << endl;

    // declare variable to store temp gradient
    // Obtain the u_l values
    double values[4];
    comp0->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    double vel_value_x = values[0];
    double vel_x_gradient_x = values[1];
    double vel_x_gradient_y = values[2];
    double vel_x_gradient_z = values[3];
    // cout << "vel_value_x: " << vel_value_x << " vel_x_gradient_x: " << vel_x_gradient_x << " vel_x_gradient_y: " << vel_x_gradient_y << " vel_x_gradient_z: " << vel_x_gradient_z << endl;

    comp1->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    double vel_value_y = values[0];
    double vel_y_gradient_x = values[1];
    double vel_y_gradient_y = values[2];
    double vel_y_gradient_z = values[3];
    // cout << "vel_value_y: " << vel_value_y << " vel_y_gradient_x: " << vel_y_gradient_x << " vel_y_gradient_y: " << vel_y_gradient_y << " vel_y_gradient_z: " << vel_y_gradient_z << endl;

    comp2->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    double vel_value_z = values[0];
    double vel_z_gradient_x = values[1];
    double vel_z_gradient_y = values[2];
    double vel_z_gradient_z = values[3];
    // cout << "vel_value_z: " << vel_value_z << " vel_z_gradient_x: " << vel_z_gradient_x << " vel_z_gradient_y: " << vel_z_gradient_y << " vel_z_gradient_z: " << vel_z_gradient_z << endl;

    if (index == 5)
      exit(0);

    // Get the values
    double *fluid_velocity_x = fluid_comp0->GetValues();
    double *fluid_velocity_y = fluid_comp1->GetValues();
    double *fluid_velocity_z = fluid_comp2->GetValues();

    double diameter = m_diameter_values[internal_level] * 1e-4; // Convert to meters
    double tau = (m_fluid_rho * diameter * diameter) / (18 * m_fluid_viscosity);
    cout << "m_particle_rho: " << m_fluid_rho << " Diameter: " << diameter << " m_fluid_viscosity: " << m_fluid_viscosity << " tau: " << tau << endl;

    double gamma = 0.001;

    // Calculate the drift velocity in x-direction
    double drift_vel_x_old = solution[index];
    double drift_velocity_x = -1.0 * (vel_value_x * vel_x_gradient_x + vel_value_y * vel_y_gradient_x + vel_value_z * vel_z_gradient_x);
    drift_velocity_x += -1.0 * (1.0 / tau) * (solution[index] - fluid_velocity_x[index]) + (1.0 - gamma) * m_g_array[0];
    cout << "Solution x: " << solution[index] << " Fluid Velocity x: " << fluid_velocity_x[index] << " gamma: " << gamma << " g: " << m_g_array[0] << " tau: " << tau << endl;
    drift_velocity_x = (drift_velocity_x * timestep) + drift_vel_x_old;

    // drift_velocity_x = fluid_velocity_x[index]; // For now, set the drift velocity to the fluid velocity.
    cout << "Fluid Velocity x: " << fluid_velocity_x[index] << " Drift Velocity x: " << drift_velocity_x << endl;

    // Calculate the drift velocity in y-direction
    double drift_vel_y_old = solution[m_n_velocity_points + index];
    double drift_velocity_y = -1.0 * (vel_value_x * vel_x_gradient_y + vel_value_y * vel_y_gradient_y + vel_value_z * vel_z_gradient_y);
    drift_velocity_y += -1.0 * (1.0 / tau) * (solution[m_n_velocity_points + index] - fluid_velocity_y[index]) + (1.0 - gamma) * m_g_array[1];
    cout << "Drift Velocity y: " << drift_velocity_y << endl;
    drift_velocity_y = (drift_velocity_y * timestep) + drift_vel_y_old;

    // drift_velocity_y = fluid_velocity_y[index]; // For now, set the drift velocity to the fluid velocity.
    cout << "Fluid Velocity y: " << fluid_velocity_y[index] << " Drift Velocity y: " << drift_velocity_y << endl;

    // Calculate the drift velocity in z-direction
    double drift_vel_z_old = solution[2 * m_n_velocity_points + index];
    double drift_velocity_z = -1.0 * (vel_value_x * vel_x_gradient_z + vel_value_y * vel_y_gradient_z + vel_value_z * vel_z_gradient_z);
    drift_velocity_z += -1.0 * (1.0 / tau) * (solution[2 * m_n_velocity_points + index] - fluid_velocity_z[index]) + (1.0 - gamma) * m_g_array[2];
    cout << "Drift Velocity z: " << drift_velocity_z << endl;
    drift_velocity_z = (drift_velocity_z * timestep) + drift_vel_z_old;

    // drift_velocity_z = fluid_velocity_z[index]; // For now, set the drift velocity to the fluid velocity.
    cout << "Fluid Velocity z: " << fluid_velocity_z[index] << " Drift Velocity z: " << drift_velocity_z << endl;

    continue;

    // Assign the computed values to the solution array
    solution[index] = drift_velocity_x;
    solution[m_n_velocity_points + index] = drift_velocity_y;
    solution[2 * m_n_velocity_points + index] = drift_velocity_z;

    // Obtain the FESpace for the drift velocity
    TFESpace3D *fespace_drift = drift_velocity_fevect->GetFESpace3D();

    // Get the Active Bound
    int N_Active = fespace_drift->GetActiveBound();
    // get Total DOF
    int N_DOF = fespace_drift->GetN_DegreesOfFreedom();

    // Number of Diriichlet DOF
    int N_DirichletDof = N_DOF - N_Active;

    // Now each array has 3 * of N_DOF values, each section for each component of velocity.
    // In one section of N_DOF, there will be N_Active + ( N_DOF - N_Active) values, where the last N_DOF - N_Active values are Dirichlet values.
    // We will set these values as zeros for each component of the drift velocity.
    memset(solution + N_Active, 0, N_DirichletDof * SizeOfDouble);
    memset(solution + 1 * N_DOF + N_Active, 0, N_DirichletDof * SizeOfDouble);
    memset(solution + 2 * N_DOF + N_Active, 0, N_DirichletDof * SizeOfDouble);

    // Free the memory allocated for the FEFunction3D
    delete comp0;
    delete comp1;
    delete comp2;

    delete fluid_comp0;
    delete fluid_comp1;
    delete fluid_comp2;
  }
}

void TSystemPBE3D::SolveDriftVelocity(double timestep, int internal_level, TFEVectFunct3D *drift_fevect, TFEVectFunct3D *particle_fevect, TFEVectFunct3D *fluid_fevect, int RK_Method)
{
  // Assume the incoming velocity is u_l ( which goes into PBE assembly)
  // We need to Solve for v(drift): which is u (fluid) + u_l (particle)
  // After solving for v(drift), we need to update the particle velocity u_l = v(drift) - u(fluid)
  // The equation to solve is d(v)/dt = -(v.grad)v - (1/\tau)(v - u) + (1 - mu)g
  // Obtain the drift velocity for the internal level
  TFEVectFunct3D *drift_velocity_fevect = drift_fevect;
  // obtain particle velocity array
  double* particle_velocity = particle_fevect->GetValues();
  double* fluid_velocity    = fluid_fevect->GetValues();
  double* drift_velocity    = drift_velocity_fevect->GetValues();


  // Store the Components of the drift velocity and the gradients
  TFEFunction3D *comp0 = drift_velocity_fevect->GetComponent(0);
  TFEFunction3D *comp1 = drift_velocity_fevect->GetComponent(1);
  TFEFunction3D *comp2 = drift_velocity_fevect->GetComponent(2);

  double *drift_velocity_x_new = new double[m_n_velocity_points]();
  double *drift_velocity_y_new = new double[m_n_velocity_points]();
  double *drift_velocity_z_new = new double[m_n_velocity_points]();

  // For ith Internal point, loop over all the physical points
  for (int index = 0; index < m_n_velocity_points; index++)
  {
    // obtain the gradients of the drift velocity at current point
    double x_coord = m_physical_coordinates[index];
    double y_coord = m_physical_coordinates[m_n_velocity_points + index];
    double z_coord = m_physical_coordinates[2 * m_n_velocity_points + index];
    int cell_id = (int)m_physical_coordinates[3 * m_n_velocity_points + index];
    TBaseCell *cell = m_fevect_drift->GetFESpace3D()->GetCollection()->GetCell(cell_id);

    // declare variable to store temp gradient
    double values[4];
    comp0->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    double vel_value_x = values[0];
    double vel_x_gradient_x = values[1];
    double vel_x_gradient_y = values[2];
    double vel_x_gradient_z = values[3];

    comp1->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    double vel_value_y = values[0];
    double vel_y_gradient_x = values[1];
    double vel_y_gradient_y = values[2];
    double vel_y_gradient_z = values[3];

    comp2->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    double vel_value_z = values[0];
    double vel_z_gradient_x = values[1];
    double vel_z_gradient_y = values[2];
    double vel_z_gradient_z = values[3];

    // Get the values
    double fluid_velocity_x = fluid_velocity[index];
    double fluid_velocity_y = fluid_velocity[m_n_velocity_points + index];
    double fluid_velocity_z = fluid_velocity[2 * m_n_velocity_points + index];

    double drift_velocity_x = drift_velocity[index];
    double drift_velocity_y = drift_velocity[m_n_velocity_points + index];
    double drift_velocity_z = drift_velocity[2 * m_n_velocity_points + index];

    m_fluid_rho = 1000; // Here fluid refers to the fluid particles dispersed in the air // THIVIN-HARDCODED
    double diameter = m_diameter_values[internal_level] * 1e-4; // Convert to meters
    double tau = (m_fluid_rho * diameter * diameter) / (18 * m_fluid_viscosity);

    // Stage 1
    double k1_x = -1.0 * (drift_velocity_x * vel_x_gradient_x + drift_velocity_y * vel_y_gradient_x + drift_velocity_z * vel_z_gradient_x);
    k1_x += -1.0 * (1.0 / tau) * (drift_velocity_x - fluid_velocity_x) + (1.0 - gamma) * m_g_array[0];

    double k1_y = -1.0 * (drift_velocity_x * vel_x_gradient_y + drift_velocity_y * vel_y_gradient_y + drift_velocity_z * vel_z_gradient_y);
    k1_y += -1.0 * (1.0 / tau) * (drift_velocity_y - fluid_velocity_y) + (1.0 - gamma) * m_g_array[1];

    double k1_z = -1.0 * (drift_velocity_x * vel_x_gradient_z + drift_velocity_y * vel_y_gradient_z + drift_velocity_z * vel_z_gradient_z);
    k1_z += -1.0 * (1.0 / tau) * (drift_velocity_z - fluid_velocity_z) + (1.0 - gamma) * m_g_array[2];

    if (RK_Method == 1)
    {
        // Update drift velocities
        drift_velocity_x_new[index] = drift_velocity_x + timestep * k1_x;
        drift_velocity_y_new[index] = drift_velocity_y + timestep * k1_y;
        drift_velocity_z_new[index] = drift_velocity_z + timestep * k1_z;

      // break out of the loop
      continue;
    }
    // Intermediate values for drift velocity components v_x, v_y, v_z
    double drift_velocity_x_mid = drift_velocity_x + 0.5 * timestep * k1_x;
    double drift_velocity_y_mid = drift_velocity_y + 0.5 * timestep * k1_y;
    double drift_velocity_z_mid = drift_velocity_z + 0.5 * timestep * k1_z;

    // Stage 2
    // Substitute all the v_x, v_y, v_z values with the intermediate values
    double k2_x = -1.0 * (drift_velocity_x_mid * vel_x_gradient_x + drift_velocity_y_mid * vel_y_gradient_x + drift_velocity_z_mid * vel_z_gradient_x);
    k2_x += -1.0 * (1.0 / tau) * (drift_velocity_x_mid - fluid_velocity_x) + (1.0 - gamma) * m_g_array[0];

    double k2_y = -1.0 * (drift_velocity_x_mid * vel_x_gradient_y + drift_velocity_y_mid * vel_y_gradient_y + drift_velocity_z_mid * vel_z_gradient_y);
    k2_y += -1.0 * (1.0 / tau) * (drift_velocity_y_mid - fluid_velocity_y) + (1.0 - gamma) * m_g_array[1];

    double k2_z = -1.0 * (drift_velocity_x_mid * vel_x_gradient_z + drift_velocity_y_mid * vel_y_gradient_z + drift_velocity_z_mid * vel_z_gradient_z);
    k2_z += -1.0 * (1.0 / tau) * (drift_velocity_z_mid - fluid_velocity_z) + (1.0 - gamma) * m_g_array[2];

    // Update drift velocities
    drift_velocity_x_new[index] = drift_velocity_x + timestep * k2_x;
    drift_velocity_y_new[index] = drift_velocity_y + timestep * k2_y;
    drift_velocity_z_new[index] = drift_velocity_z + timestep * k2_z;
  }

  // Obtain the FESpace for the drift velocity
  TFESpace3D *fespace_drift = drift_velocity_fevect->GetFESpace3D();

  // Get the Active Bound
  int N_Active = fespace_drift->GetActiveBound();
  // get Total DOF
  int N_DOF = fespace_drift->GetN_DegreesOfFreedom();

  // Number of Diriichlet DOF
  int N_DirichletDof = N_DOF - N_Active;

  // memcpy the new drift velocity values to the solution array
  memcpy(drift_velocity, drift_velocity_x_new, m_n_velocity_points * SizeOfDouble);
  memcpy(drift_velocity + m_n_velocity_points, drift_velocity_y_new, m_n_velocity_points * SizeOfDouble);
  memcpy(drift_velocity + 2 * m_n_velocity_points, drift_velocity_z_new, m_n_velocity_points * SizeOfDouble);


  // Compute the L2 norm of the difference between the old and new drift velocity
  double l2_norm_drift_velocity_x = L2Norm(m_n_velocity_points, drift_velocity);
  double l2_norm_drift_velocity_y = L2Norm(m_n_velocity_points, drift_velocity + m_n_velocity_points);
  double l2_norm_drift_velocity_z = L2Norm(m_n_velocity_points, drift_velocity + 2 * m_n_velocity_points);

  // cout << "L2 Norm of Drift Velocity x: " << l2_norm_drift_velocity_x << " Drift Velocity y: " << l2_norm_drift_velocity_y << " Drift Velocity z: " << l2_norm_drift_velocity_z << endl;

  // Particle Velocity (u_l) -> Drift Velocity (v) - Fluid Velocity (u)
  // Update the particle velocity
  for (int i = 0; i < m_n_velocity_points; i++)
  {
    particle_velocity[i]                            = fluid_velocity[i];
    particle_velocity[m_n_velocity_points + i]      = fluid_velocity[m_n_velocity_points + i];
    particle_velocity[2 * m_n_velocity_points + i]  = fluid_velocity[2 * m_n_velocity_points + i];
  }

  // Also set the particle velocity also to be zero on the boundaries, so what when we add them ,t he boundary values becomes zero
  memset(particle_velocity + N_Active, 0, N_DirichletDof * SizeOfDouble);
  memset(particle_velocity + 1 * N_DOF + N_Active, 0, N_DirichletDof * SizeOfDouble);
  memset(particle_velocity + 2 * N_DOF + N_Active, 0, N_DirichletDof * SizeOfDouble);

  // Now each array has 3 * of N_DOF values, each section for each component of velocity.
  // In one section of N_DOF, there will be N_Active + ( N_DOF - N_Active) values, where the last N_DOF - N_Active values are Dirichlet values.
  // We will set these values as zeros for each component of the drift velocity.
  memset(drift_velocity + N_Active, 0, N_DirichletDof * SizeOfDouble);
  memset(drift_velocity + 1 * N_DOF + N_Active, 0, N_DirichletDof * SizeOfDouble);
  memset(drift_velocity + 2 * N_DOF + N_Active, 0, N_DirichletDof * SizeOfDouble);

  // Free the memory allocated for the FEFunction3D
  delete comp0;
  delete comp1;
  delete comp2;

  // Free the memory created for the new drift velocity values
  delete[] drift_velocity_x_new;
  delete[] drift_velocity_y_new;
  delete[] drift_velocity_z_new;
}

// Helper functions
double ComputeConvectiveTerm(
    int index,
    const double *u1, const double *u2, const double *u3,            // Velocity components
    const double *grad_x, const double *grad_y, const double *grad_z // Gradients
)
{
  return u1[index] * grad_x[index] +
         u2[index] * grad_y[index] +
         u3[index] * grad_z[index];
}

// Calculate RHS of equation for one component
double ComputeComponentRHS(
    int index, int component,
    const double *u_x, const double *u_y, const double *u_z,
    const double *V_x, const double *V_y, const double *V_z,
    const double *grad_u_x, const double *grad_u_y, const double *grad_u_z,
    const double *grad_V_x, const double *grad_V_y, const double *grad_V_z,
    double tau, double gamma, const double *g_array)
{
  // (u·∇)u term
  double conv_u = ComputeConvectiveTerm(index,
                                        u_x, u_y, u_z,
                                        grad_u_x, grad_u_y, grad_u_z);
  // (u·∇)V term
  double conv_uV = ComputeConvectiveTerm(index,
                                         u_x, u_y, u_z,
                                         grad_V_x, grad_V_y, grad_V_z);

  // (V·∇)u term
  double conv_Vu = ComputeConvectiveTerm(index,
                                         V_x, V_y, V_z,
                                         grad_u_x, grad_u_y, grad_u_z);

  // (V·∇)V term
  double conv_V = ComputeConvectiveTerm(index,
                                        V_x, V_y, V_z,
                                        grad_V_x, grad_V_y, grad_V_z);

  // Get appropriate component of V and g
  const double *V_comp = (component == 0) ? V_x : (component == 1) ? V_y : V_z;

  return -conv_u - conv_uV - conv_Vu - conv_V - (1.0 / tau) * V_comp[index] + (1.0 - gamma) * g_array[component];

}

void TSystemPBE3D::SolveDriftVelocity(double timestep, int internal_level,
                                      TFEVectFunct3D *fluid_fevect,    // Fluid velocity FE vector function
                                      TFEVectFunct3D *particle_fevect, // Drift velocity FE vector function
                                      double *particle_velocity, double *fluid_velocity)
{
  /* Solving the component equation:
   * ∂u/∂t + ∂V/∂t + (u·∇)u + (u·∇)V + (V·∇)u + (V·∇)V = -1/τ(V) + (1-γ)g
    Here "V" is the u_l (particle velocity) and "u" is the fluid velocity.
    Here ∂u/∂t = 0, as the fluid velocity is constant.
   */

  // Pre-allocate arrays for all gradients and values
  double *u_x = new double[m_n_velocity_points]();
  double *u_x_dx = new double[m_n_velocity_points]();
  double *u_x_dy = new double[m_n_velocity_points]();
  double *u_x_dz = new double[m_n_velocity_points]();

  double *u_y = new double[m_n_velocity_points]();
  double *u_y_dx = new double[m_n_velocity_points]();
  double *u_y_dy = new double[m_n_velocity_points]();
  double *u_y_dz = new double[m_n_velocity_points]();

  double *u_z = new double[m_n_velocity_points]();
  double *u_z_dx = new double[m_n_velocity_points]();
  double *u_z_dy = new double[m_n_velocity_points]();
  double *u_z_dz = new double[m_n_velocity_points]();

  double *V_x = new double[m_n_velocity_points]();
  double *V_x_dx = new double[m_n_velocity_points]();
  double *V_x_dy = new double[m_n_velocity_points]();
  double *V_x_dz = new double[m_n_velocity_points]();

  double *V_y = new double[m_n_velocity_points]();
  double *V_y_dx = new double[m_n_velocity_points]();
  double *V_y_dy = new double[m_n_velocity_points]();
  double *V_y_dz = new double[m_n_velocity_points]();

  double *V_z = new double[m_n_velocity_points]();
  double *V_z_dx = new double[m_n_velocity_points]();
  double *V_z_dy = new double[m_n_velocity_points]();
  double *V_z_dz = new double[m_n_velocity_points]();

  // Get all components
  TFEFunction3D *fluid_comp0 = fluid_fevect->GetComponent(0);
  TFEFunction3D *fluid_comp1 = fluid_fevect->GetComponent(1);
  TFEFunction3D *fluid_comp2 = fluid_fevect->GetComponent(2);

  TFEFunction3D *particle_velocity_comp0 = particle_fevect->GetComponent(0);
  TFEFunction3D *particle_velocity_comp1 = particle_fevect->GetComponent(1);
  TFEFunction3D *particle_velocity_comp2 = particle_fevect->GetComponent(2);

  // Pre-calculate all gradients
  for (int index = 0; index < m_n_velocity_points; index++)
  {
    double x_coord = m_physical_coordinates[index];
    double y_coord = m_physical_coordinates[m_n_velocity_points + index];
    double z_coord = m_physical_coordinates[2 * m_n_velocity_points + index];
    int cell_id = (int)m_physical_coordinates[3 * m_n_velocity_points + index];
    TBaseCell *cell = particle_fevect->GetFESpace3D()->GetCollection()->GetCell(cell_id);

    double values[4];

    // Get all fluid velocity gradients
    fluid_comp0->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    u_x[index] = values[0];
    u_x_dx[index] = values[1];
    u_x_dy[index] = values[2];
    u_x_dz[index] = values[3];

    fluid_comp1->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    u_y[index] = values[0];
    u_y_dx[index] = values[1];
    u_y_dy[index] = values[2];
    u_y_dz[index] = values[3];

    fluid_comp2->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    u_z[index] = values[0];
    u_z_dx[index] = values[1];
    u_z_dy[index] = values[2];
    u_z_dz[index] = values[3];

    // Get all drift velocity gradients
    particle_velocity_comp0->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    V_x[index] = values[0];
    V_x_dx[index] = values[1];
    V_x_dy[index] = values[2];
    V_x_dz[index] = values[3];

    particle_velocity_comp1->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    V_y[index] = values[0];
    V_y_dx[index] = values[1];
    V_y_dy[index] = values[2];
    V_y_dz[index] = values[3];

    particle_velocity_comp2->FindGradientLocal(cell, cell_id, x_coord, y_coord, z_coord, values);
    V_z[index] = values[0];
    V_z_dx[index] = values[1];
    V_z_dy[index] = values[2];
    V_z_dz[index] = values[3];

  }

  // Arrays for new drift velocity components
  double *particle_velocity_x_new = new double[m_n_velocity_points]();
  double *particle_velocity_y_new = new double[m_n_velocity_points]();
  double *particle_velocity_z_new = new double[m_n_velocity_points]();

  // Calculate physical parameters
  double diameter = m_diameter_values[internal_level] * 1e-4; // Convert to meters
  double tau = (m_particle_rho[internal_level] * diameter * diameter) / (18 * m_fluid_viscosity);
  double gamma = m_fluid_rho / m_particle_rho[internal_level];

  // Main computation loo
  // Main computation loop - First Order Euler
  for (int index = 0; index < m_n_velocity_points; index++)
  {
      // Calculate RHS for all components
      double dV_x_dt = ComputeComponentRHS(
          index,                    // Current point index
          0,                        // X component indicator
          u_x, u_y, u_z,           // Fluid velocity components
          V_x, V_y, V_z,           // Current drift velocities
          u_x_dx, u_x_dy, u_x_dz,  // Gradients of fluid velocity x-component
          V_x_dx, V_x_dy, V_x_dz,  // Gradients of drift velocity x-component
          tau,                      // Relaxation time
          gamma,                    // Density ratio
          m_g_array                 // Gravity vector
      );

      double dV_y_dt = ComputeComponentRHS(
          index,
          1,                        // Y component indicator
          u_x, u_y, u_z,
          V_x, V_y, V_z,
          u_y_dx, u_y_dy, u_y_dz,  // Gradients of fluid velocity y-component
          V_y_dx, V_y_dy, V_y_dz,  // Gradients of drift velocity y-component
          tau,
          gamma,
          m_g_array
      );

      double dV_z_dt = ComputeComponentRHS(
          index,
          2,                        // Z component indicator
          u_x, u_y, u_z,
          V_x, V_y, V_z,
          u_z_dx, u_z_dy, u_z_dz,  // Gradients of fluid velocity z-component
          V_z_dx, V_z_dy, V_z_dz,  // Gradients of drift velocity z-component
          tau,
          gamma,
          m_g_array
      );

      // Simple Euler update
      particle_velocity_x_new[index] = V_x[index] + timestep * dV_x_dt;
      particle_velocity_y_new[index] = V_y[index] + timestep * dV_y_dt;
      particle_velocity_z_new[index] = V_z[index] + timestep * dV_z_dt;

  }
  // Update solution arrays
  memcpy(particle_velocity, particle_velocity_x_new, m_n_velocity_points * SizeOfDouble);
  memcpy(particle_velocity + m_n_velocity_points, particle_velocity_y_new, m_n_velocity_points * SizeOfDouble);
  memcpy(particle_velocity + 2 * m_n_velocity_points, particle_velocity_z_new, m_n_velocity_points * SizeOfDouble);

  // // Get the Active DOF from particle velocity FESpace
  // TFESpace3D* fespace_particle = particle_fevect->GetFESpace3D();
  // int N_Active = fespace_particle->GetActiveBound();
  // int N_DOF = fespace_particle->GetN_DegreesOfFreedom();

  // // Set the Dirichlet values to zero
  // memset(particle_velocity + N_Active, 0, (N_DOF - N_Active) * SizeOfDouble);
  // memset(particle_velocity + m_n_velocity_points + N_Active, 0, (N_DOF - N_Active) * SizeOfDouble);
  // memset(particle_velocity + 2 * m_n_velocity_points + N_Active, 0, (N_DOF - N_Active) * SizeOfDouble);

  // Compute L-inf norm of the difference between old and new particle velocity
  double l2_norm_particle_velocity_x = LinfNorm(m_n_velocity_points, particle_velocity);

  // Loop over particle_velocity and find the minimum value
  double min_particle_velocity_x = particle_velocity[0];
  for (int i = 1; i < m_n_velocity_points; i++)
  {
    if (particle_velocity[i] < min_particle_velocity_x)
    {
      min_particle_velocity_x = particle_velocity[i];
    }
  }

  cout << "-----------------------------------------------------------------------------------------" << endl;
  cout << "---- L2 Norm of Particle Velocity x: " << l2_norm_particle_velocity_x << endl;
  cout << "---- Min Particle Velocity x: " << min_particle_velocity_x << endl;
  cout << "-----------------------------------------------------------------------------------------" << endl;

  // Clean up
  delete[] particle_velocity_x_new;
  delete[] particle_velocity_y_new;
  delete[] particle_velocity_z_new;

  // Clean up gradient arrays
  delete[] u_x;
  delete[] u_x_dx;
  delete[] u_x_dy;
  delete[] u_x_dz;
  delete[] u_y;
  delete[] u_y_dx;
  delete[] u_y_dy;
  delete[] u_y_dz;
  delete[] u_z;
  delete[] u_z_dx;
  delete[] u_z_dy;
  delete[] u_z_dz;
  delete[] V_x;
  delete[] V_x_dx;
  delete[] V_x_dy;
  delete[] V_x_dz;
  delete[] V_y;
  delete[] V_y_dx;
  delete[] V_y_dy;
  delete[] V_y_dz;
  delete[] V_z;
  delete[] V_z_dx;
  delete[] V_z_dy;
  delete[] V_z_dz;

  // Clean up FE components
  delete fluid_comp0;
  delete fluid_comp1;
  delete fluid_comp2;
  delete particle_velocity_comp0;
  delete particle_velocity_comp1;
  delete particle_velocity_comp2;
}

void TSystemPBE3D::AssembleMRhs()
{
  // this is set to true for direct solver factorization
  factorize = true;

  int i, N_DOF_low, N_Active;

  for (i = Start_Level; i < N_Levels; i++)
  {
    N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
    N_Active = FeSpaces[i]->GetActiveBound();

    /** initialize matrices and rhs */
    MMatRhsAssemble[i]->Reset();

    // assemble
    MMatRhsAssemble[i]->Assemble3D();

    /** free the Mass mat array, no need in time loop */
    MMatRhsAssemble[i]->DeAllocate();

    /** set rhs for Dirichlet nodes */
    memcpy(SolArray[i] + N_Active, RhsArray[i] + N_Active, (N_DOF_low - N_Active) * SizeOfDouble);
  } //  for(i=Start_Level;i<N_Levels;i++)

} // TSystemMatScalar3D::AssembleMRhs

void TSystemPBE3D::AssembleARhs()
{
  // this is set to true for direct solver factorization
  factorize = true;

  int i, N_DOF_low, N_Active;

  for (i = Start_Level; i < N_Levels; i++)
  {
    N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
    N_Active = FeSpaces[i]->GetActiveBound();

    /** reset the matrix and rhs */
    AMatRhsAssemble[i]->Reset();

    // assemble
    AMatRhsAssemble[i]->Assemble3D();

    /** set rhs for Dirichlet nodes */
    memcpy(SolArray[i] + N_Active, RhsArray[i] + N_Active, (N_DOF_low - N_Active) * SizeOfDouble);
  } //   for(i=Start_Level;i<N_Le

} // TSystemMatScalar3D::AssembleARhs

void TSystemPBE3D::AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol
#ifdef _MPI
                                   ,
                                   double **Rhs_array
#endif
)
{
  int i, N_Active;
  double tau;

  if (SystMatAssembled)
  {
    OutPut("System is has to be restored before calling AssembleSystMat! " << endl);
    exit(0);
  }

  SQMATRICES[0] = sqmatrixM[N_Levels - 1];

  N_Active = FeSpaces[N_Levels - 1]->GetActiveBound();
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  memset(B, 0, N_DOF * SizeOfDouble);

  /** old rhs multiplied with current subtime step and theta3 on B */
  Daxpy(N_Active, tau * TDatabase::TimeDB->THETA3, oldrhs, B);

  /** add rhs from current sub time step to rhs array B */
  Daxpy(N_Active, tau * TDatabase::TimeDB->THETA4, rhs, B);

  /** M = M + (- tau*THETA2)A */
  MatAdd(sqmatrixM[N_Levels - 1], sqmatrixA[N_Levels - 1], -tau * TDatabase::TimeDB->THETA2);
  gamma = -tau * TDatabase::TimeDB->THETA2; // set current factor of steady state matrix

  /** defect = M * oldsol */
  memset(defect, 0, N_DOF * SizeOfDouble);
  MatVectActive(sqmatrixM[N_Levels - 1], oldsol, defect);
  // cout << "defect " << Ddot(N_Active, sol, sol)<< endl;

  /** B:= B + defec  */
  Daxpy(N_Active, 1, defect, B);

  /** set Dirichlet values */
  memcpy(B + N_Active, rhs + N_Active, (N_DOF - N_Active) * SizeOfDouble);
  memcpy(sol + N_Active, rhs + N_Active, (N_DOF - N_Active) * SizeOfDouble);

  /** assemble the system matrix */
  for (i = Start_Level; i < N_Levels; i++)
  {
    if (i == N_Levels - 1)
    {
      MatAdd(sqmatrixM[i], sqmatrixA[i], -gamma + tau * TDatabase::TimeDB->THETA1);
    }
    else
    {
      MatAdd(sqmatrixM[i], sqmatrixA[i], tau * TDatabase::TimeDB->THETA1);
    }

#ifdef _MPI
    SQMATRICES[0] = sqmatrixM[i];
#endif
  }
  gamma = tau * TDatabase::TimeDB->THETA1;

// have to shift this in pardirectsolver
#ifdef _OMPONLY
  if (SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
    DS->AssembleMatrix(sqmatrixM[N_Levels - 1]);
#endif

  SystMatAssembled = TRUE;

} // AssembleSystMat

void TSystemPBE3D::RestoreMassMat()
{
  int i;

  if (SystMatAssembled)
  {
    // restore the mass matrix
    for (i = Start_Level; i < N_Levels; i++)
      MatAdd(sqmatrixM[i], sqmatrixA[i], -gamma);

    gamma = 0.;
    SystMatAssembled = FALSE;
  }
  else
  {
    cout << "System is not assembled to restore " << endl;
    exit(0);
  }
}

void TSystemPBE3D::Solve(double *sol)
{
  switch (SOLVER)
  {
  case AMG_SOLVE:
    Solver(sqmatrixM[N_Levels - 1], B, sol);
    break;

  case GMG:
    if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
    {
      memcpy(Itmethod_sol, sol, N_DOF * SizeOfDouble);
      memcpy(Itmethod_rhs, B, N_DOF * SizeOfDouble);
    }
    else
    {
      Itmethod_sol = sol;
      Itmethod_rhs = B;
    }

    /** solve linear system */
    Itmethod->Iterate(sqmatrices, NULL, Itmethod_sol, Itmethod_rhs);
#ifdef _MPI
    if (TDatabase::ParamDB->SC_SMOOTHER_SCALAR == 6)
      ParComm[N_Levels - 1]->CommUpdateH2(Itmethod_sol);
#endif
    if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
    {
      memcpy(sol, Itmethod_sol, N_DOF * SizeOfDouble);
    }
    break;

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
#endif
    // this is set to false for direct solver factorization
    factorize = false;
    break;

  default:
    OutPut("Unknown Solver" << endl);
    exit(4711);
    ;
  }
}

void TSystemPBE3D::Solve_Pardiso(double *sol, int iter_num)
{
  switch (SOLVER)
  {

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
    // DirectSolver(sqmatrixM[N_Levels-1], B, sol);
    PardisoDirectSolverWithObject(sqmatrixM[N_Levels - 1], B, sol, iter_num, pardiso_solver);
#endif
    // this is set to false for direct solver factorization
    factorize = false;
    break;

  default:
    OutPut("Unknown Solver" << endl);
    exit(4711);
    ;
  }
}

double TSystemPBE3D::GetResidual(double *sol)
{
  double residual_scalar = 0.0;

  if (SystMatAssembled)
  {
    memset(defect, 0, N_DOF * SizeOfDouble);
    ScalarDefect(sqmatrixM[N_Levels - 1], sol, B, defect, residual_scalar);

#ifdef _MPI
    residual_scalar = 0.0;
    double sum = 0.;
    int i, rank;
    MPI_Comm_rank(Comm, &rank);
    int *master = ParComm[N_Levels - 1]->GetMaster();
    for (i = 0; i < N_DOF; i++)
    {
      if (master[i] != rank)
        continue;
      residual_scalar += defect[i] * defect[i];
    }
    MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
    residual_scalar = sqrt(sum);
#endif
  }
  else
  {
    OutPut("Assemble the System Matrix before calculating the GetResidual" << endl);
    exit(4711);
    ;
  }
  return residual_scalar;
}
