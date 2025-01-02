/** ************************************************************************ 
*
* @class     TSystemPBE3D
* @brief     stores the information of a timedependent part of a 3D Population Balance
* @author    Thivin Anandh
* @date      08-oct-2024
* @History 
 ************************************************************************  */


#ifndef __SYSTEMPBE3D__
#define __SYSTEMPBE3D__

#include <SquareMatrix3D.h>
#include <SystemCD3D.h>

/**class for 3D scalar system matrix */
class TSystemPBE3D : public TSystemCD3D
{
  protected:
#ifdef _MPI
    TParDirectSolver *TDS;
#endif   
    /** M mass matrix */
    TSquareMatrix3D **sqmatrixM;
    
    /** working rhs, used in AssembleSystMat() */
    double *B;   
   
    /** to store defect */
    double *defect;   
    
    /** factor that multplied with Mat A in working rhs */
    double gamma;   
    
    bool factorize;

    
    /** instance of the Assemble class */
    TAssembleMat3D **MMatRhsAssemble;
    
//     /** Stiffness part of the SUPG matrix */
//     TSquareMatrix3D *sqmatrixK;    
//     
//     /** time-consistent part of the SUPG matrix */
//     TSquareMatrix3D *sqmatrixS;
    
    /** Discrete form of the M and rhs matrics */
    TDiscreteForm3D *DiscreteFormMRhs, *DiscreteFormRhs; 
    
    /** Systmat assemble indicator */
    bool SystMatAssembled;
    
  public:
    /** constructor */
     TSystemPBE3D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver);

    /** destrcutor */
    ~TSystemPBE3D();

    /** methods */
    void Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue, TAuxParam3D *aux);

    /** Function, which calls the Assembly function with the NSE3D Values */
    void Init_with_NSEValues(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue, TAuxParam3D *aux);

    /** Function, initialises the parameters for Drift velocity */
    void Init_for_drift_velocity(TFEVectFunct3D** fevect_drift_array, double* g, double fluid_rho, double* particle_rho, double fluid_viscosity,
          int N_internal, double* diameter_values, TFEVectFunct3D* fevect_b, int n_velocity_points);
    
    /** Solve for Drift Velocity */
    void SolveDriftVelocity(double timestep, int i, double *sol);
    
    // Values
    TFEVectFunct3D **m_fevect_drift_array; // Stores the drift velocity array
    double* m_g_array;
    double m_fluid_rho;  // Stores fluid density
    double* m_particle_rho; // Stores particle density for each internal point
    double m_fluid_viscosity; // Stores fluid viscosity
    int m_N_internal; // Stores the number of internal points
    double* m_internal_values; // Stores the internal values (Number Concentration of particles)
    double* m_diameter_values; // Stores the diameter values of the particles
    TFEVectFunct3D* m_fevect_fluid_velocity;
    int m_n_velocity_points;

    double* m_physical_coordinates; // to store all physical co-ordinates where the fluid velocity is stored
    TFEVectFunct3D* m_fevect_drift; // Fefunction to access the physical co-ordinates where the fluid velocity is stored

    /** return the stiffness matric */
    TSquareMatrix3D **GetAMatrix()
    { return sqmatrixA; }
    
    /** assemble the Mass mat and rhs */
    void AssembleMRhs(); 
    
    /** assemble the stifness mat and rhs */
    void AssembleARhs();   
    
//     /** M = M + (tau*THETA1)*A */ 
//     /** B = (tau*THETA1)*rhs +(tau*THETA2)*oldrhs + [ M - (tau*THETA2)A]*oldsol */  
    void AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol
#ifdef _MPI
                         , double **Rhs_array
#endif
                         );
    
    /** restoring the mass matrix */
    void RestoreMassMat();
    
     /** solve the system matrix */
    void Solve(double *sol);  
    
    /** return the residual of the system for the given sol*/
    double GetResidual(double *sol);
    
    double value(double *sol,int N){
      int i;
      double sum=0.;
      for(i=0;i<N;i++)	sum+=sol[i];
      return sum;
    }
    
};

#endif

