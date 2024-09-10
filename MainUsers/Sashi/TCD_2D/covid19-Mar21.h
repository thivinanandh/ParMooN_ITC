
void ExampleFile()
{
 #ifdef _MPI
 if(TDatabase::ParamDB->Par_P0==1)
 #endif 
  OutPut("Example: Covid19-example.h" << endl) ;
 
  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;
  
  #define __LD__
  // #define __LV__  
  // #define __LA__  

  #define __LOCKDOWNMODEL__
  #define __OSNODALPT__
}

// initial conditon
void InitialValues(int N_Coord, double *X, double *values)
{
 double k=1;
 double t=TDatabase::TimeDB->CURRENTTIME ;
 double x = X[0];
 double y = X[1];    
 double la = X[2];  
 double lv = X[3];
 double InitialPopulation = X[4];
 double ld = X[5]; 
 double vsymp = TDatabase::ParamDB->P1; 
 double f5, lvsim2 = (lv-vsymp)*(lv-vsymp);
 double k5 = 2.3937;
 double dist, a1= 4.71, b1= 0.28, c1 = 0.12, b2= 1.62, a2 = 0.42 ;

 if(lv<vsymp)
  { f5 = k5*exp( -lvsim2/(2*(vsymp/3.)*(vsymp/3.)) );  }
 else
  { f5 = k5*exp( -lvsim2/(2*((1.-vsymp)/3.)*((1.-vsymp)/3.)) ); }
 
 if(fabs(ld)<1e-4)
  {
   values[0] = 0.;
  }
 else
  {
   dist = (a1*exp(- (la-b1)*(la-b1)/(c1*c1)))*(1./(ld*a2*sqrt(2*Pi)))*exp(-pow( (log(ld)-b2), 2.0)/(2*a2*a2))*f5;    
   values[0] =  InitialPopulation*dist;
  }
 //gnuplot plot [200:350] (1./(x*0.42*sqrt(2*22/7)))*exp(-(log(x)-1.61)* (log(x)-1.61))

  // values[0] = 60;

  //  if(TDatabase::ParamDB->P2==1)
  // if(fabs(dist)>1)
  //  cout<< la << " : " << lv <<" : " << ld <<   " InitialPopulation  : " << InitialPopulation << endl;
}


// ====================================================================================
// functions for L0 system
// ====================================================================================
void BoundCondition_LminLMax_L0(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;

  // cond_Lmin = DIRICHLET;
  // cond_Lmax = DIRICHLET;

}

void Set_sigma_SD(double t, double &NoNucTime, double offsetsize, double delta)
{
  double sd=TDatabase::ParamDB->P13;
  double bsd = TDatabase::ParamDB->P14;  //0.1
  double sd_c = TDatabase::ParamDB->P15;  //0.5
  double f1;

 //  set Interactive index
//  if(t<15) // from 23rd Mar to 7th April
//    {
//     f1 = .125;
//    }
//    else if (t<21) //till 15th April
//    {
//     f1 = .125 - (t - 15)*0.025/7;
//    }
//    else if(t<36) // till 30th April
//    {
//     f1 = 0.1 - (t - 21)*0.005/15;
//    }
//    else if(t<51) // till 15th May
//    {
//     f1 = 0.095 - (t - 36)*0.03/15;
//    }   
//    else
//    {
//     f1 = 0.065;
//    }  

  f1 = 1.;

  // if(fabs(t- NoNucTime)<=0.5 ) // sunday break means, no nuc from Sat 8pm to Mon 6 am!
  // { 
  //  f1 = 0.25;
  // }

  TDatabase::ParamDB->P7 = f1;

 //  set Social Distancing  index

  // Mar 2020 fit
  // if(t<15)  
  //  { sd = 0.7 + t*0.02/15.; }  // from 23rd Mar to 7th April
  // else if(t<36.) 
  //  { sd = 0.72 + (t-15.)*0.09/21.; }   // till 30th April 
  // else if(t<51.) 
  //  { sd = 0.81 - (t-36.)*0.06/15.; }  // till 15th May
  // else if(t<72.)
  //  { sd = 0.75 + (t-51.)*0.0125/21; }  // till  
  // else if(t<102.)
  //  {      
  //   sd = 0.7625 + (t-72.)*0.025/30; 
  //  }
  // else if(t<132.)
  //  {      
  //   sd = 0.7875 - (t-102.)*0.0125/30; 
  //  }
  // else  
  //  {      
  //   sd = 0.775 + (t - 132.)*delta/offsetsize; 
  //  }

  // sd = 0.727; 
  //Mar21: 15 days lockdown
  // if(t<35)  
  //   {  sd = 0.727; }  // from 23rd Mar to 26th Apr,  // no lockdown
  // else if(t<50) 
  //  { sd = 0.727 + (t-35.)*0.055/15.; }   // till   
  // else if(t<65.) 
  //  { sd = 0.782 -(t-50.)*0.007/15.;  }
  // else
  //  { sd = 0.775; }

  //Mar21: 21 days lockdown
  // if(t<35)  
  //   {  sd = 0.727; }  // from 23rd Mar to 26th Apr,  // no lockdown
  // else if(t<56) 
  //  { sd = 0.727 + (t-35.)*0.0775/21.; }   // till   
  // else if(t<71.) 
  //  { sd = 0.8045 -(t-56.)*0.0295/15.;  } //unlock
  // else
  //  { sd = 0.775; }

  //Mar21: 30 days lockdown
  if(t<35)  
    {  sd = 0.727; }  // from 23rd Mar to 26th Apr,  // no lockdown
  else if(t<65) 
   { sd = 0.727 + (t-35.)*0.116/30.; }   // till   
  else if(t<80.) 
   { sd = 0.843 -(t-65.)*0.068/15.;  } //unlock
  else
   { sd = 0.775; }

   // Sep 2020 fit
      // if(t<30.)
      //  { sd = 0.77 + t*0.002/15; }
      // else if (t<153)
      //  { sd = 0.774 - (t-30)*0.003/15; } //till end of 2020
      // else
      //  { sd = 0.749 + (t-153)*delta/offsetsize; } 


 TDatabase::ParamDB->P13 = 1. -  1./(1.+ exp(-(sd - sd_c)/bsd)); // f3
  // 1. -  1./(1.+ exp(-( (0.3 + 0.125*sin(2.*22.*t/(30.*7.)) + (t-75.)*0.035/15.)   - 0.5)/0.1))
  //  OutPut("f1: " << TDatabase::ParamDB->P7 << " " );
  //  OutPut("f3: " << TDatabase::ParamDB->P13<<endl);
}

void GrowthAndB_Nuc_L0(int N_Inputs, double *Inn, double *Out) 
{
  double B_Nuc;
  double R_0  = TDatabase::ParamDB->P6; // Nucleation factor
  double f1  = TDatabase::ParamDB->P7; // interactive index
  // double bsigma  = TDatabase::ParamDB->P8;     
  // double sigma_c = TDatabase::ParamDB->P9;  
  double vsymp = TDatabase::ParamDB->P1;
  // double health  = TDatabase::ParamDB->P10;    
  // double bhealth = TDatabase::ParamDB->P11;  
  // double health_c = TDatabase::ParamDB->P12;  
  double f3  = TDatabase::ParamDB->P13; 
  double a4= 4.7, b4=0.28, c4_sqr=0.12*0.12;      
  // double t=TDatabase::TimeDB->CURRENTTIME ; 

  if(N_Inputs!=6)
  {
    printf("N_Inputs != 5,  Get_GrowthAndB_Nuc: %d \n",N_Inputs); 
#ifdef _MPI
     MPI_Finalize();
#endif       
    exit(0);        
  }

  // double x = Inn[0];  
  //  y = Inn[1];  
  double la = Inn[2]; // age
  double lv = Inn[3]; // 
  double IntL_val = Inn[4];
  double SuscepPopRatio = Inn[5];

  double k5 = 2.3937;



  // if(t<20)
  // { sd = 0.1 + t*0.2/20; } // f3 =0.89
  // else if(t<75)
  // { sd = 0.3; }  
  // else if(t<150)
  // { sd = 0.3 + (t-75.)*0.2/75; } // f3 =0.5  
  // else 
  // { sd = 0.5; } // f3 =0.5


       // plot [0:0.5] 1. -  1./(1.+ exp(-(x - 0.5)/0.1)) //for gnuplot

  //  f1 = 1./(1.+ exp(-(sigma - sigma_c)/bsigma));
  // f1 = 1.;
  // f2 = 1. -  1./(1.+ exp(-(health - health_c)/bhealth));
  // f2 = 1.;

 
  double f5, f4 = a4*exp(-(la -b4)* (la -b4)/c4_sqr);

  double lvsim2 = (lv-vsymp)*(lv-vsymp);
  if(lv<vsymp)
   {
    f5 = k5*exp( -lvsim2/(2.*(vsymp/3.)*(vsymp/3.)) );
   }
   else
   {
    f5 = k5*exp( -lvsim2/(2*((1.-vsymp)/3.)*((1.-vsymp)/3.)) );   
   }
  // f5 = 1;

  // R = R_0*f1*f2*f3*f4*f5;

  B_Nuc = SuscepPopRatio*R_0*f1*f3*f4*f5*IntL_val;  // 
  //  B_Nuc =  R_0*f1*f3*IntL_val;  // 
  //  if(B_Nuc<0)
  //  {
    // cout << "B_Nuc: " << B_Nuc << " IntL_val "  <<IntL_val << " : " << R <<endl;
  //   // exit(0);
  //  }  

  Out[0] = 1; // growth is 1
  Out[1] = B_Nuc; 
  Out[2] = 0; //L_Max Bound value,  if DIRICHLET 
    // cout << "B_Nuc: " <<  Out[1] <<endl;
}

void BilinearCoeffs_L0(int n_points, int N_Dim, double **Coords,
                             double **parameters, double **coeffs)
{
  int i;
  double eps_L, *coeff;
  // double inn[5], out[3];
  double ld;
  double t=TDatabase::TimeDB->CURRENTTIME ;
  double k=1;
  double b_nuc=0., ekt = exp(-k*t);
  // inn[0] = Coords[0][0]; //x
  // inn[1] = Coords[1][0]; //y
  // inn[2] = Coords[2][0]; // la
  // inn[3] = Coords[3][0]; // lv
  // inn[4] = Coords[4][0]; // lv
  // int Cell_NO = (int)Coords[4][1];

  if(TDatabase::ParamDB->REACTOR_P3)
    eps_L = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps_L = 0.;

  //  if(Cell_NO == 0)
  //  {
  //    GrowthAndB_Nuc_L0(5, inn, out) ;
  //    b_nuc = out[1];
  //   //  cout << "BilinearCoeffs_L0 B_nuc: " << b_nuc <<endl;
  //  }

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // ld = Coords[5][i]; 

    // diffusion
    coeff[0] = eps_L;
    
    // convection in l1 direction (-ve ==> growth term will be the convection)
    coeff[1] = -1;
    
    // reaction term
    coeff[2] = 0;
    
    //rhs  
    // coeff[3] =  b_nuc;
    coeff[3] = 0.;
    // coeff[3] =  -k*(exp(-k*t))*sin(Pi*ld)*cos(Pi*lv)*cos(Pi*la) 
    //              + Pi*(exp(-k*t))*cos(Pi*ld)*cos(Pi*lv)*cos(Pi*la);
    // coeff[3] =  -k*ekt*sin(Pi*ld)  + Pi*ekt*cos(Pi*ld) + eps_L*Pi*Pi*sin(Pi*ld);
    // coeff[3] =  -k*ekt*sin(Pi*ld)  + eps_L*Pi*Pi*sin(Pi*ld);
        // cout << "ld: " << ld << " "<< coeff[3]  <<endl;
  }
}

// ====================================================================================
// functions for L1 system
// ====================================================================================
 void BoundCondition_LminLMax_L1(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
  // cond_Lmin = DIRICHLET;
  // cond_Lmax = DIRICHLET;  
}


void GrowthAndB_Nuc_L1(int N_Inputs, double *Inn, double *Out) 
{

  double la, G;
  double K_g = TDatabase::ParamDB->REACTOR_P20;  //G_DimLessFact  
  double arisk  = TDatabase::ParamDB->P5; // Age_risk
  double p  = 2; // Age_offsetPower

  if(N_Inputs!=5)
  {
    printf("N_Inputs != 5,  Get_GrowthAndB_Nuc: %d \n",N_Inputs); 
#ifdef _MPI
     MPI_Finalize();
#endif       
    exit(0);        
  }
  
  la = Inn[2]; // age
  G =  K_g*pow((la - arisk), p);

  Out[0] = G; // growth
  Out[1] = 0;  //B_Nuc no nucleation
  Out[2] = 0; //L_Max Bound value,  if DIRICHLET 
  //  cout <<"Lv growth: " <<Out[0]<<endl;
}

void BilinearCoeffs_L1(int n_points, int N_Dim, double **Coords,
                             double **parameters, double **coeffs)
{
  int i;
  double eps_L, *coeff;
 
  // if(TDatabase::ParamDB->REACTOR_P3)
  //   eps_L = 1.0/TDatabase::ParamDB->REACTOR_P3;
  // else
    eps_L = 0.;


  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    // diffusion
    coeff[0] = 0.;
    
    // convection in l1 direction (-ve ==> growth term will be the convection)
    coeff[1] = -1;
    
    // reaction term
    coeff[2] = 0;
    
    //rhs  
    coeff[3] =  0;
  }
}
// ====================================================================================
// functions for L3 system
// ====================================================================================
 
void BoundCondition_LminLMax_L2(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
//   cond_Lmin = DIRICHLET;
//   cond_Lmax = DIRICHLET;  
}

void GrowthAndB_Nuc_L2(int N_Inputs, double *Inn, double *Out) 
{
  Out[0] = 0.; // no growth,  L_Min Bound value, if DIRICHLET 
  Out[1] = 0.; // no nucleation
  Out[2] = 0; //L_Max Bound value,  if DIRICHLET   
}

void BilinearCoeffs_L2(int n_points, int N_Dim, double **Coords,
                             double **parameters, double **coeffs)
{
  int i;
  double eps_L, *coeff;
 
  // if(TDatabase::ParamDB->REACTOR_P3)
  //   eps_L = 1.0/TDatabase::ParamDB->REACTOR_P3;
  // else
    eps_L = 0.;


  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    // diffusion
    coeff[0] = 0;
    
    // convection in l1 direction (-ve ==> growth term will be the convection)
    coeff[1] = -1;
    
    // reaction term
    coeff[2] = 0;
    
    //rhs  
    coeff[3] =  0;
  }
}

void Gamma_Q(int N_Coord, double *X, double *values)
{
  double bv = TDatabase::ParamDB->P0;
  double vsymp = TDatabase::ParamDB->P1;
  double bd = TDatabase::ParamDB->P2;
  double dsymp = TDatabase::ParamDB->P3;
  double ba = TDatabase::ParamDB->P4;
  double arisk = TDatabase::ParamDB->P5;
  
//  double t=TDatabase::TimeDB->CURRENTTIME ;
//  double x = X[0];
//  double y = X[1];    
 double la = X[2];  
 double lv = X[3];
 double ld = X[5];  

  //  cout<< bv << " vsymp  : " <<vsymp<<": " <<  bd << " dsymp  : " <<dsymp <<": "<<  ba << " values  : " <<arisk << endl;

//  values[0] = (1./(1.+ exp(-2.*bv*(lv-vsymp))))*(1./(1.+ exp(-2.*bd*(ld-dsymp))))*(1./(1.+ exp(-2.*ba*(la-arisk))));
//  values[0] =  0.75*(1./(1.+ exp(-(ld-dsymp)/2.))) ;

//  double val1 = 1./(1.+ exp(-(ld-dsymp)/bd)); // qurantine
//  double fact1 =   2.*(t/75.) + 30.*(1.-t/75.); 
//  double val2 = (1./(1.+ exp(-(ld-dsymp)/fact1))); // incubation period, the spread will be more initally

//  double theta = t/(30.+ld); // testing kid also be less compare to infected population over time
//  double fact = theta*2.0 + (1.-theta)*0.01; // initilly, the nucleation will be more and over a period people become moree H and SD
  
  if(ld<1)
  { values[0] = 1.; }
  else
   { 
    // values[0] = 0.9*(1./(1.+ exp(-(ld-dsymp)/bd)));  
    values[0] = 0.9*(1./(1.+ exp(-(lv-vsymp)/bv)))*(1./(1.+ exp(-(ld-dsymp)/bd)));
   }
//  values[0] = 0.6;

  //  if(values[0]<0.)
  //  cout<< ld << " values  : " << values[0] << endl;
}

void GetReactionFactor(int *N_LnodalPos, double **LnodalPos, double *ReactionFactor)
{
  int i,j,k;
  double *ValL0, *ValL1L0;
  double ld, lv, la, C_R, r_la, r_lv, f1_R, C_ID, f1_Id, f2_Id;

  // double bv = TDatabase::ParamDB->P0;
  // double vsymp = TDatabase::ParamDB->P1;
  // double bd = TDatabase::ParamDB->P2;
  // double dsymp = TDatabase::ParamDB->P3;  
  // double arisk  = TDatabase::ParamDB->P5; // Age_risk
  // double kr = 0.3;
  // double kid = 0.05;
  int N_L1L0 =  N_LnodalPos[0]*N_LnodalPos[1];
  double t=TDatabase::TimeDB->CURRENTTIME ; 

  for(i=0;i<N_LnodalPos[2]; ++i)
   {
    ValL1L0 = ReactionFactor + 2*i*N_L1L0;
    la = LnodalPos[2][i];
     
     // r_la = (la - arisk)*(la - arisk);
    for(j=0;j<N_LnodalPos[1]; ++j)
     {
      ValL0 = ValL1L0 + 2*j*N_LnodalPos[0];
      lv = LnodalPos[1][j];

      //  r_lv = 1./(1.+ exp(-2.*bv*(lv-vsymp)));
      for(k=0;k<N_LnodalPos[0]; ++k)
      {
       ld = LnodalPos[0][k];

       if(lv<0.64)
        { 
         f1_Id = 0;
        }
       else
        {
      // Mar 2020 fit          
      //    if(t<21)
      //    {  
      //     f1_Id = 0.0475 - t*0.0075/21.;
      //    }
      //   else if(t<36.) 
      //    {
      //     f1_Id = 0.0475 + (t-21)*0.0075/36.;

      //    }
      //   else if(t<90) 
      //   {
      //    f1_Id = 0.05 - (t-36.)*0.03/55.; 
      //   }
      //  else if(t<160) 
      //   {
      //    f1_Id = 0.02 - (t-90.)*0.003/25.;
      //   }
      //  else 
      //   {
      //    f1_Id = 0.008;
      //   }

      // // Sep 2020 fit
      // if(t<30) 
      //   {
      //    f1_Id = 0.02 - t*0.012/30.;
      //   }
      //  else 
      //   {
      //    f1_Id = 0.007;
      //   }  

         f1_Id = 0.007;

       } //  if(lv>=0.64)
        
      // Mar 2020 fit
      //   if(t<65)  //  till 26th May
      //   {
      //    f1_R = 0.01 + t*0.0375/65.;
      //   }   
      //  else if(t<120) 
      //   {
      //   //  f1_R = 0.0475; // 17th June model
      //   f1_R = 0.0475 + (t-65.)*0.005/60; // 3 Aug model
      //   }      
      //  else  
      //   {
      //    f1_R = 0.0475; // 17th June model
      //   f1_R = 0.0521 + 0.005/75; // 15 Aug model
      //   f1_R = 0.06 + t*0.00215/15; // 5th Sep model 1
      //   }  

         // Sep 2020 model
         if(t<30)
          {  f1_R =  0.04 + 1./(1+exp(-0.03*t))*0.75/15; }
         else
          { f1_R =  0.076; } 
 

       //Mar 2021
       if(t<30)
         { 
          f1_R =   0.04 + 1./(1+exp(-0.03*t))*0.45/15;
         }
       else
        {
         f1_R = 0.076;
        }   

       C_R = f1_R;
       C_ID = f1_Id;

       ValL0[2*k] = C_R+C_ID;
       ValL0[2*k+1] = f1_R;       
        // cout<< " ld: " <<ValL0[k] <<endl;
      }
     }
   }
}

 void GetLExampleData(int *L_NVertices, double *L_StartEnd, BoundCond1D **BoundConLminLMax, 
                      DoubleFunctND **GrowthAndB_Nuc)
  {
   // L_d  
   L_StartEnd[0] = TDatabase::ParamDB->REACTOR_P12;
   L_StartEnd[1] = TDatabase::ParamDB->REACTOR_P13;
   L_NVertices[0] = (int)TDatabase::ParamDB->REACTOR_P11;
   GrowthAndB_Nuc[0] = GrowthAndB_Nuc_L0;
   BoundConLminLMax[0] = BoundCondition_LminLMax_L0;

   // L_v  
   L_StartEnd[2] = 0.;
   L_StartEnd[3] = 1.;
   L_NVertices[1] = (int)TDatabase::ParamDB->REACTOR_P15;
   GrowthAndB_Nuc[1] = GrowthAndB_Nuc_L1;
   BoundConLminLMax[1] = BoundCondition_LminLMax_L1;
 
   // L_a  
   L_StartEnd[4] = 0.;
   L_StartEnd[5] = 1.0; // la/125
   L_NVertices[2] = (int)TDatabase::ParamDB->REACTOR_P16;
   GrowthAndB_Nuc[2] = GrowthAndB_Nuc_L2;
   BoundConLminLMax[2] = BoundCondition_LminLMax_L2;
  }

// ====================================================================================
// functions for spatial system
// ====================================================================================

void BoundCondition_Spatial(int BD_ID, double t, BoundCond &cond)
{
  switch(BD_ID)
  {
   case 0:
   case 1:
   case 2:
   case 3:
     cond = NEUMANN;
   break;

   default: 
      cout << "wrong boundary part number" << endl;
      exit(0);
   break;
  }
}

void BoundValue_Spatial(int BdComp, double Param, double &value)
{
  value =0;
}

void BilinearCoeffs_Spatial(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double *coeff, *param;
  double x, y;
 
  double eps;
  if(TDatabase::ParamDB->REACTOR_P2)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P2;
  else
    eps = 0.;
 
  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
     }
    else
     {
      coeff[1] = 0.;  // u1
      coeff[2] = 0.;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0; // f 
  }
}
