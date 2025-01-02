/** ************************************************************************ 
* @brief     source file for TSystemPBE3D
* @author    Sashikumaar Ganesan
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

TSystemPBE3D::TSystemPBE3D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver)
                          :TSystemCD3D(N_levels, fespaces, sol, rhs, disctype, solver)
{
  int i;
  
  /** M mass and system matrix */
  sqmatrixM = new TSquareMatrix3D*[N_Levels];
  N_Matrices++;
  
  MMatRhsAssemble = new TAssembleMat3D*[N_Levels];
  
  for(i=Start_Level;i<N_Levels;i++)
   {
    sqmatrixM[i] = new TSquareMatrix3D(sqstructure[i]);       
   }
  
  /** working rhs, used in AssembleSystMat() */
  B = new double[N_DOF];
  defect = new double[N_DOF];
  
  gamma =0.;
  
//   /** time-consistent part of the SUPG matrix */
//   if(Disctype==SDFEM || Disctype==SUPG)
//    {
//     sqmatrixS = new TSquareMatrix3D(sqstructure); 
//     N_Matrices++;
//     
//     sqmatrixK = new TSquareMatrix3D(sqstructure); 
//     N_Matrices++; 
//    }
   
  SystMatAssembled  = FALSE;
    
} // constructor


TSystemPBE3D::~TSystemPBE3D()
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


void TSystemPBE3D::Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue,
                              TAuxParam3D *aux )
{
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
  BoundaryValues[0] = BoundValue;
  
  TDiscreteForm3D *DiscreteFormMRhs_Galerkin;
  TDiscreteForm3D *DiscreteFormARhs_Galerkin; 
//   TDiscreteForm3D *DiscreteFormMRhs_SUPG;
//   TDiscreteForm3D *DiscreteFormARhs_SUPG;

  
  InitializeDiscreteFormsScalar(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormRhs, BilinearCoeffs);
  
    switch(Disctype)
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
} // Init


// For Population Balance Equation with NSE Values. 
// IT alsi includes the C . \grad(u_p)
void TSystemPBE3D::Init_with_NSEValues(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue,
                              TAuxParam3D *aux )
{
  cout << "Correct Function for Init_with_NSEValues" << endl;
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
  BoundaryValues[0] = BoundValue;
  
  TDiscreteForm3D *DiscreteFormMRhs_Galerkin;
  TDiscreteForm3D *DiscreteFormARhs_Galerkin; 
//   TDiscreteForm3D *DiscreteFormMRhs_SUPG;
//   TDiscreteForm3D *DiscreteFormARhs_SUPG;

  
  InitializeDiscreteFormsScalarNSE(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormRhs, BilinearCoeffs);
  
    switch(Disctype)
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
} // Init


// Function to import values for the drift velocity
void TSystemPBE3D::Init_for_drift_velocity(TFEVectFunct3D** fevect_drift_array, double* g, double fluid_rho, double* particle_rho, double fluid_viscosity,
          int N_internal, double* diameter_values, TFEVectFunct3D* fevect_b, int n_velocity_points)
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

  // Set up Array and memory to store the gradients of the drift velocity
  double** particle_velocity_gradient_x = new double*[m_N_internal];
  double** particle_velocity_gradient_y = new double*[m_N_internal];
  double** particle_velocity_gradient_z = new double*[m_N_internal];

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
  m_physical_coordinates = new double [4 * m_n_velocity_points]();  // here the 4th dimension will be used to store cell id.
  TFESpace3D* fespace_drift = m_fevect_drift->GetFESpace3D();

  // Create a new FeVectFunction to Store the Physical Co-ordinates
  TFEVectFunct3D* fevect_physical_coordinates = new TFEVectFunct3D(fespace_drift, "Physical_Coordinates", "Physical_Coordinates", m_physical_coordinates, m_n_velocity_points, 3);
  fevect_physical_coordinates->GridToDataWithCellid(); // populates cell id on 4th dimension


} 

void TSystemPBE3D::SolveDriftVelocity(double timestep, int internal_level, double* solution)
{
  // Assign the solution to the solution array
  TFEVectFunct3D* drift_velocity_fevect = m_fevect_drift_array[internal_level];

  // For ith Internal point, loop over all the physical points
  for (int index = 0 ; index < m_n_velocity_points ; index++)
  {
    // obtain the gradients of the drift velocity at current point
    double x_coord = m_physical_coordinates[index];
    double y_coord = m_physical_coordinates[m_n_velocity_points + index];
    double z_coord = m_physical_coordinates[2 * m_n_velocity_points + index];
    int cell_id = (int)m_physical_coordinates[3 * m_n_velocity_points + index];
    TBaseCell* cell = m_fevect_drift->GetFESpace3D()->GetCollection()->GetCell(cell_id);
    
    // Store the Components of the drift velocity and the gradients
    TFEFunction3D* comp0 = drift_velocity_fevect->GetComponent(0);
    TFEFunction3D* comp1 = drift_velocity_fevect->GetComponent(1);
    TFEFunction3D* comp2 = drift_velocity_fevect->GetComponent(2);

    
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

    // Store the components
    TFEFunction3D* fluid_comp0 = m_fevect_fluid_velocity->GetComponent(0);
    TFEFunction3D* fluid_comp1 = m_fevect_fluid_velocity->GetComponent(1);
    TFEFunction3D* fluid_comp2 = m_fevect_fluid_velocity->GetComponent(2);

    // Get the values
    double* fluid_velocity_x = fluid_comp0->GetValues();
    double* fluid_velocity_y = fluid_comp1->GetValues();
    double* fluid_velocity_z = fluid_comp2->GetValues();

    // Check if there is a mismatch on the interpolated values and the actual values
    // This should not be occurring, since we are picking the solution values at nodal points, and we are interpolating at the nodal points. 
    // if (( abs(solution[index] - vel_value_x) > 1e-5  ) || (abs(solution[m_n_velocity_points + index] - vel_value_y) > 1e-5) || (abs(solution[2 * m_n_velocity_points + index] - vel_value_z) > 1e-5))
    // {
    //   cout << "Mismatch in the values of the drift velocity at index : " << index << " Internal Level : " << internal_level << endl;
    //   cout << " Solution Value x: " << solution[index] << " Interpolated Value : " << vel_value_x ;
    //   cout << " Solution Value y: " << solution[m_n_velocity_points + index] << " Interpolated Value : " << vel_value_y ;
    //   cout << " Solution Value z: " << solution[2 * m_n_velocity_points + index] << " Interpolated Value : " << vel_value_z ;

    //   vel_value_x = 0;
    //   vel_value_y = 0;
    //   vel_value_z = 0;

    //   vel_x_gradient_x = 0;
    //   vel_x_gradient_y = 0;
    //   vel_x_gradient_z = 0;

    //   vel_y_gradient_x = 0;
    //   vel_y_gradient_y = 0;
    //   vel_y_gradient_z = 0;
    // }

    double diameter = m_diameter_values[internal_level] * 1e-6; // Convert to meters
    double tau = (m_particle_rho[internal_level] * diameter * diameter) / (18 * m_fluid_viscosity);

    double gamma = 0.001;

    // Calculate the drift velocity in x-direction
    double drift_velocity_x = -1.0 * (vel_value_x* vel_x_gradient_x + vel_value_y * vel_y_gradient_x + vel_value_z * vel_z_gradient_x) ;
    drift_velocity_x += -1.0 * (1.0/tau) * (solution[index] - fluid_velocity_x[index]) + (1.0 - gamma) * m_g_array[0];
    drift_velocity_x = drift_velocity_x * timestep + vel_value_x;
    
    drift_velocity_x = fluid_velocity_x[index]; // For now, set the drift velocity to the fluid velocity.

    // Calculate the drift velocity in y-direction
    double drift_velocity_y = -1.0 * (vel_value_x* vel_x_gradient_y + vel_value_y * vel_y_gradient_y + vel_value_z * vel_z_gradient_y) ;
    drift_velocity_y += -1.0 * (1.0/tau) * (solution[m_n_velocity_points + index] - fluid_velocity_y[index]) + (1.0 - gamma) * m_g_array[1];
    drift_velocity_y = drift_velocity_y * timestep + vel_value_y;

    drift_velocity_y = fluid_velocity_y[index]; // For now, set the drift velocity to the fluid velocity.

    // Calculate the drift velocity in z-direction
    double drift_velocity_z = -1.0 * (vel_value_x* vel_x_gradient_z + vel_value_y * vel_y_gradient_z + vel_value_z * vel_z_gradient_z) ;
    drift_velocity_z += -1.0 * (1.0/tau) * (solution[2 * m_n_velocity_points + index] - fluid_velocity_z[index]) + (1.0 - gamma) * m_g_array[2];
    drift_velocity_z = drift_velocity_z * timestep + vel_value_z;

    drift_velocity_z = fluid_velocity_z[index]; // For now, set the drift velocity to the fluid velocity.

    // Assign the computed values to the solution array
    solution[index] = drift_velocity_x;
    solution[m_n_velocity_points + index] = drift_velocity_y;
    solution[2 * m_n_velocity_points + index] = drift_velocity_z;

    // Free the memory allocated for the FEFunction3D
    delete comp0;
    delete comp1;
    delete comp2;

    delete fluid_comp0;
    delete fluid_comp1;
    delete fluid_comp2;

  }
}




void TSystemPBE3D::AssembleMRhs()
{
  //this is set to true for direct solver factorization
  factorize = true;
  
  int i, N_DOF_low, N_Active;

  for(i=Start_Level;i<N_Levels;i++)
  {      
      N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
      N_Active =  FeSpaces[i]->GetActiveBound();
      
      /** initialize matrices and rhs */
      MMatRhsAssemble[i]->Reset(); 
      
      // assemble
      MMatRhsAssemble[i]->Assemble3D();
      
      /** free the Mass mat array, no need in time loop */
      MMatRhsAssemble[i]->DeAllocate();
      
      
      /** set rhs for Dirichlet nodes */
      memcpy(SolArray[i]+N_Active, RhsArray[i]+N_Active, (N_DOF_low - N_Active)*SizeOfDouble);         
    } //  for(i=Start_Level;i<N_Levels;i++)

} // TSystemMatScalar3D::AssembleMRhs 


void TSystemPBE3D::AssembleARhs()
{
  //this is set to true for direct solver factorization
  factorize = true;
  
  int i, N_DOF_low, N_Active;

   for(i=Start_Level;i<N_Levels;i++)
    {    
     N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
     N_Active =  FeSpaces[i]->GetActiveBound();

     /** reset the matrix and rhs */
     AMatRhsAssemble[i]->Reset(); 
    
     // assemble
     AMatRhsAssemble[i]->Assemble3D();     

     /** set rhs for Dirichlet nodes */
     memcpy(SolArray[i]+N_Active, RhsArray[i]+N_Active, (N_DOF_low - N_Active)*SizeOfDouble);           
    }//   for(i=Start_Level;i<N_Le  
    
} // TSystemMatScalar3D::AssembleARhs 

void TSystemPBE3D::AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol
#ifdef _MPI
                                             , double **Rhs_array
#endif
                                             )
{
    int i, N_Active;
    double tau;
    
    if(SystMatAssembled)
     {
      OutPut("System is has to be restored before calling AssembleSystMat! " <<endl);
      exit(0);
     }
    
    SQMATRICES[0] = sqmatrixM[N_Levels-1];

    N_Active =  FeSpaces[N_Levels-1]->GetActiveBound();     
    tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;  
    
    memset(B, 0, N_DOF*SizeOfDouble); 

    /** old rhs multiplied with current subtime step and theta3 on B */
    Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3,  oldrhs, B);    

    /** add rhs from current sub time step to rhs array B */
    Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4,  rhs, B);    

    /** M = M + (- tau*THETA2)A */
     MatAdd(sqmatrixM[N_Levels-1], sqmatrixA[N_Levels-1], - tau*TDatabase::TimeDB->THETA2);
     gamma = -tau*TDatabase::TimeDB->THETA2;  // set current factor of steady state matrix

     /** defect = M * oldsol */
     memset(defect, 0, N_DOF*SizeOfDouble);  
     MatVectActive(sqmatrixM[N_Levels-1], oldsol, defect); 
    //cout << "defect " << Ddot(N_Active, sol, sol)<< endl;      
    
     /** B:= B + defec  */
     Daxpy(N_Active, 1, defect, B);

     /** set Dirichlet values */
     memcpy(B+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);  
     memcpy(sol+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
          
     /** assemble the system matrix */
     for(i=Start_Level;i<N_Levels;i++)   
      {
       if(i==N_Levels-1)
         { MatAdd(sqmatrixM[i], sqmatrixA[i], -gamma + tau*TDatabase::TimeDB->THETA1);}
        else
         { MatAdd(sqmatrixM[i], sqmatrixA[i], tau*TDatabase::TimeDB->THETA1);} 
         
#ifdef _MPI  
       SQMATRICES[0] = sqmatrixM[i];  
#endif
       }
     gamma = tau*TDatabase::TimeDB->THETA1;
     
//have to shift this in pardirectsolver     
#ifdef _OMPONLY     
    if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
      DS->AssembleMatrix(sqmatrixM[N_Levels-1]);
#endif
     
     SystMatAssembled  = TRUE;

} // AssembleSystMat

void TSystemPBE3D::RestoreMassMat()
{
 int i;

  if(SystMatAssembled)
   {
     // restore the mass matrix
     for(i=Start_Level;i<N_Levels;i++)  
      MatAdd(sqmatrixM[i], sqmatrixA[i], -gamma);
     
     gamma = 0.;
     SystMatAssembled  = FALSE;
   }
  else
  {
    cout << "System is not assembled to restore " <<endl;
    exit(0);
  }

}

void TSystemPBE3D::Solve(double *sol)
{  
    switch(SOLVER)
     {
      case AMG_SOLVE:
         Solver(sqmatrixM[N_Levels-1], B, sol);
      break;

      case GMG:
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          memcpy(Itmethod_sol, sol, N_DOF*SizeOfDouble);
          memcpy(Itmethod_rhs, B, N_DOF*SizeOfDouble);
         }
        else
         {
          Itmethod_sol = sol;
          Itmethod_rhs = B;
         }
      
         /** solve linear system */
        Itmethod->Iterate(sqmatrices, NULL, Itmethod_sol, Itmethod_rhs);
#ifdef _MPI
    if(TDatabase::ParamDB->SC_SMOOTHER_SCALAR==6)
         ParComm[N_Levels-1]->CommUpdateH2(Itmethod_sol);
#endif
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          memcpy(sol, Itmethod_sol, N_DOF*SizeOfDouble);
         }
      break;

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
#endif
	//this is set to false for direct solver factorization
        factorize = false;
      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }
     
}

double TSystemPBE3D::GetResidual(double *sol)
{
  double residual_scalar=0.0;
  
  if(SystMatAssembled)
   {
    memset(defect, 0, N_DOF*SizeOfDouble);         
    ScalarDefect(sqmatrixM[N_Levels-1], sol, B, defect, residual_scalar);
    
#ifdef _MPI 
    residual_scalar = 0.0;
    double sum =0.;
    int i,rank;
    MPI_Comm_rank(Comm, &rank);
    int *master = ParComm[N_Levels-1]->GetMaster();
    for(i=0;i<N_DOF;i++)
    {
      if(master[i]!=rank)    continue;
      residual_scalar += defect[i]*defect[i];
    }
   MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
   residual_scalar = sqrt(sum);
#endif
   }
  else
   {
    OutPut("Assemble the System Matrix before calculating the GetResidual" << endl);
    exit(4711);;   
   }
   return residual_scalar;    
}



