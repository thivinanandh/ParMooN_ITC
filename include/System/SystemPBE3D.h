/** ************************************************************************ 
*
* @class     TSystemPBE3D
* @brief     stores the information of a scalar 3D+3DSurface  PDE system 
* @author    Sashikumaar Ganesan
* @date      07.10.2019
* @History 
 ************************************************************************  */


#ifndef __SYSTEMPBE3D__
#define __SYSTEMPBE3D__

#include <SquareMatrix3D.h>
#include <SystemTCD3D.h>
#include<FESpace1D.h>
#include<FEFunction1D.h>
#include <Output3D.h>

/**class for 3D scalar system matrix */
class TSystemPBE3D : public TSystemTCD3D
{
  protected:
#ifdef _MPI
  TParDirectSolver *TDS;
#endif       

  /** Domain for the internal space */
  TDomain *Domain_Intl;
  TDomain *SurfDomain;
  TCollection *Coll_Intl, *Surf_Coll;
  int *Cell_IDX, *Joint_IDX;

  /** Physical space parameters for internal space */
  int N_SpatialPts, SpatialSpace_Level;
  double *SpatialPts, *Spatial_X, *Spatial_Y, *Spatial_Z;
  TFEVectFunct3D *SpatialPos;

  /** Internal space parameters for internal space */
  int N_SurfCells;
  double *InternalPts;
  TFEVectFunct3D *InternalPos;

  /** parameters internal space multigrid solver*/
  int LEVELS_INTL, ORDER_INTL;   
    
  /** FE spaces and functionns for internal domain*/
  TFESpace1D *FeSpaces_Intl1D;  
  TFEFunction1D *FeFunction_Intl1D;
  TFEFunction3D **SpatialFEFunctions; 
  TFESpace2D *FeSpaces_Intl;

  /** size of internal space */
  int N_DoF_Intl, N_DoF_Intl1D;
  double *Sol_Intl1D, *sol_intl, *rhs_intl, *BE_MatScale;  
    
  /** assembling the convective term in physical system*/  
  TSquareMatrix3D *SQMATRICES_RTE[14];  
  TSquareMatrix3D **B1, **B2, **B3;
  TSquareMatrix3D **S1_SUPG, **S2_SUPG, **S3_SUPG, **S1S2_SUPG, **S1S3_SUPG, **S2S3_SUPG, **M1_SUPG, **M2_SUPG, **M3_SUPG;
  double gammab0, gammab1, gammab2, gammab3, gammab4, gammab5, gammab6, gammab7, gammab8, gammab9;
  double gammab10, gammab11, gammab12;
  double **Rhs1Array,  **Rhs2Array, **Rhs3Array, *rhs1old, *rhs2old, *rhs3old;
  double **Rhs_SUPG[9];

  TOutput3D *Output;

  bool factorize;
    
  /** instance of the Assemble class */
  TAssembleMat3D **RhsAssemble;
    
//  /** Stiffness part of the SUPG matrix */
//     TSquareMatrix3D *sqmatrixK;    
//     
//     /** time-consistent part of the SUPG matrix */
//     TSquareMatrix3D *sqmatrixS;
    
    /** Discrete form of the M and rhs matrics */
    // TDiscreteForm3D *DiscreteFormMRhs, *DiscreteFormRhs; 
    
  /** Systmat assemble indicator */
  bool SystMatAssembled;

#ifdef _SMPI  
  const int root = 0;
  int rank, mpi_size;

 int N_AllSurfCell;
 int *DispArray;
 int *N_Cells_ALL; 
#endif



  public:
    /** constructor */
     TSystemPBE3D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver,
                       TFEFunction3D **spatialFEFunctions, DoubleFunctND *getKernel, DoubleFunctND *getRhs, TDomain* Domain_Intl);

    /** destrcutor */
    ~TSystemPBE3D();

    /** methods */
    void Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue, 
              BoundCondFunct3D *BoundCond_Intl, BoundValueFunct3D *BoundValue_Intl,TAuxParam3D *aux,
              BoundCondFunct2D *SurfBDCond, BoundValueFunct2D * SurfBDValue);
    
    /** return the stiffness matric */
    // TSquareMatrix3D **GetAMatrix()
    // { return sqmatrixA; }
    void GetQuadIntegral(int i, int &N_Points, double * &weights, double *x, double *y, double *z,
                                 double *n1, double *n2, double *n3);
    /** assemble the internal (Mass) mat  */
    void AssembleIntlMat(int CellNo, double tau, double *XPos, double &fact, double *InternalScaling); 
    
    /** assemble the stifness mat and rhs */
    void AssembleABRhs();   

    /** assemble the stifness rhs */    
    void AssembleRhs();

//     /** M = M + (tau*THETA1)*(B1 + B2 + B3) ) */ 
   void AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol, 
                                 double b1, double b2, double b3, CoeffFct3D *SystemRhs
#ifdef _MPI
                         , double **Rhs_array
#endif
                       );
    
   void Interpolate(DoubleFunct3D *Exact, double *Sol);


    /** Get internal Points*/
    double *GetIntlPoints()
    {  return InternalPts;  }

    /** assume tha nodal values are point vales */
    int GetN_IntlPoints()
    { return N_DoF_Intl; }

    /** restoring the mass matrix */
    void RestoreMassMat();
    
    /** solve the system matrix */
    void Solve(double *sol);  
    
    /** return the internal (scalled, vector) matrix */
    double *GetIntlMat()
     {return BE_MatScale;} 

    /** kernel funcion */
    DoubleFunctND *GetKernel, *GetRhs;

    /** Compute Error */
    void GetErrors(DoubleFunctVect *Exact, TFEFunction3D *ScalarFunction, double *sol_all, double &l2);
    
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

