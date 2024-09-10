/** ************************************************************************ 
*
* @class     TSystemCD3D_3D
* @brief     stores the information of a multidimensional part of a 3D scalar system 
* @author    Sashikumaar Ganesan
* @date      06.10.2019
* @History 
 ************************************************************************  */


#ifndef __SYSTEMTCD3D_3D__
#define __SYSTEMTCD3D_3D__

#include <SquareMatrix3D.h>
#include <SystemTCD3D.h>

/**class for 3D scalar system matrix */
class TSystemTCD3D_3D : public TSystemTCD3D
{
  protected:
#ifdef _MPI
  TParDirectSolver *TDS;
#endif   
  /** Dimension of the Internal Space */
  int InternalSpaceDim;
    
  /** Domain for the internal space */
  TDomain *Domain_Intl;
  TCollection *Coll_Intl;

  /** Physical space parameters for internal space */
  int N_PhysicalPts, PhycicalSpace_Level;
  double *PhysicalPts, *Physical_X, *Physical_Y, *Physical_Z;
  TFEVectFunct3D *PhysicalPos;

  /** Internal space parameters for internal space */
  int InternalSpace_Level, N_InternalPts;
  double *InternalPts, *Internal_X, *Internal_Y, *Internal_Z;
  TFEVectFunct3D *InternalPos;

  /** parameters internal space multigrid solver*/
  int LEVELS_INTL, mg_type_intl, mg_level_intl, ORDER_INTL;   
    
  /** FE spaces and finctionns for internal domain*/
  TFESpace3D **FeSpaces_Intl;  
  TFEFunction3D *FeFunction_Intl, **FeFunctions_Intl; 

  /** size of internal space */
  int N_DoF_Intl;
  double *sol_intl, *rhs_intl, **Sol_array_Intl, **Rhs_array_Intl;  
    
    // bool factorize;

    
    /** instance of the Assemble class */
    // TAssembleMat3D **MMatRhsAssemble;
    
//     /** Stiffness part of the SUPG matrix */
//     TSquareMatrix3D *sqmatrixK;    
//     
//     /** time-consistent part of the SUPG matrix */
//     TSquareMatrix3D *sqmatrixS;
    
    /** Discrete form of the M and rhs matrics */
    // TDiscreteForm3D *DiscreteFormMRhs, *DiscreteFormRhs; 
    
    /** Systmat assemble indicator */
    // bool SystMatAssembled;
    
  public:
    /** constructor */
     TSystemTCD3D_3D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver, int intspacedim);

    /** destrcutor */
    ~TSystemTCD3D_3D();

    /** methods */
    void Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue, 
              BoundCondFunct3D *BoundCond_Intl, BoundValueFunct3D *BoundValue_Intl,TAuxParam3D *aux);
    
    /** return the stiffness matric */
    // TSquareMatrix3D **GetAMatrix()
    // { return sqmatrixA; }
    
    /** assemble the Mass mat and rhs */
    // void AssembleMRhs(); 
    
    /** assemble the stifness mat and rhs */
    // void AssembleARhs();   
    
//     /** M = M + (tau*THETA1)*A */ 
//     /** B = (tau*THETA1)*rhs +(tau*THETA2)*oldrhs + [ M - (tau*THETA2)A]*oldsol */  
    // void AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol
#ifdef _MPI
                         , double **Rhs_array
#endif
                        //  );
    
    /** restoring the mass matrix */
    // void RestoreMassMat();
    
     /** solve the system matrix */
    // void Solve(double *sol);  
    
    /** return the residual of the system for the given sol*/
    // double GetResidual(double *sol);
    
    // double value(double *sol,int N){
      // int i;
      // double sum=0.;
      // for(i=0;i<N;i++)	sum+=sol[i];
      // return sum;
    // }
    
};

#endif

