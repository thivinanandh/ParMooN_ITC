
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

  #define __AGEGROUP__

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
 double f4, f5, lvsim2 = (lv-vsymp)*(lv-vsymp);
//  double k5 = 911.58171714;
 double k5;
//  double k5 = 12.511*2.3937;
 double dist, a1= 4.71, b1= 0.28, c1 = 0.12, b2= 1.62, a2 = 0.42 ;

 
  //Jun 2021, based on data: https://ncdc.gov.in/dashboard.php
  // if(la<=10./125.) {
  //   f4 = 0.0336; }
  // else if(la<=20./125.) {
  //   f4 = 0.0841; }
  // else if(la<=30./125.){
  //   f4 = 0.218; }
  // else if(la<=40./125.){
  //   f4 = 0.2193; }    
  // else if(la<=50./125.){
  //   f4 = 0.1722; }    
  // else if(la<=60./125.){
  //   f4 = 0.1405; }    
  // else if(la<=70./125.){
  //   f4 = 0.0866; }    
  // else if(la<=80./125.){
  //   f4 = 0.0354; }    
  // else if(la<=90./125.){
  //   f4 = 0.0092; }    
  // else {
  //   f4 = 0.0012;}
  
 //needed in Int_B_Nuc()

 double norm = 1.0;
 
 f4 = 0.;

 if (TDatabase::ParamDB->PBE_P10<0.)
  {
      norm = TDatabase::ParamDB->PBE_P5*(.11/125) +
                TDatabase::ParamDB->PBE_P6*((17.-11.)/125.) + 
                TDatabase::ParamDB->PBE_P7*((44.-17.)/125.) + 
                TDatabase::ParamDB->PBE_P8*((59.-44.)/125.) + 
                TDatabase::ParamDB->PBE_P9*((125.-59.)/125.); 

  if(la<=11./125.) {
    f4 = 2.56919365672*TDatabase::ParamDB->PBE_P5/norm; 
    }
  else if(la<=17./125.) {
    f4 = 4.70419796961*TDatabase::ParamDB->PBE_P6/norm; 
    }
  else if(la<=44./125.){
    f4 =  1.04537706409*TDatabase::ParamDB->PBE_P7/norm;
     }
  else if(la<=59./125.){
    f4 = 1.88168560379*TDatabase::ParamDB->PBE_P8/norm;
     }    
  else {
    f4 = 0.42765485306*TDatabase::ParamDB->PBE_P9/norm;   
    }  
    // cout << " Norm " << norm << endl;
    k5 = 911.58171714;   
    // k5 =  886.515872968; // constant fixed to match the total numbers 
  }
 else
  { 
   f4 = 1.;
   k5 = 926.744363172;   // constant fixed to match the total numbers 
  }   

 TDatabase::ParamDB->FS_U = norm;

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
   // model used till May 2021
   //  dist = (a1*exp(- (la-b1)*(la-b1)/(c1*c1)))*(1./(ld*a2*sqrt(2*Pi)))*exp(-pow( (log(ld)-b2), 2.0)/(2*a2*a2))*f5;    

   //jun 2021, f4 - based on data: https://ncdc.gov.in/dashboard.php
   dist = f4*(1./(ld*a2*sqrt(2*Pi)))*exp(-pow( (log(ld)-b2), 2.0)/(2*a2*a2))*f5; 
   values[0] =  InitialPopulation*dist/TDatabase::ParamDB->REACTOR_P13;   
  }
  
 //gnuplot plot [200:350] (1./(x*0.42*sqrt(2*22/7)))*exp(-(log(x)-1.61)* (log(x)-1.61))

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

void Set_sigma_SD(double t, double NoNucTime, double offsetsize, double delta)
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
// 	Nationwide lockdown:
//     Phase 1: 25 March 2020 – 14 April 2020 (21 days)
//     Phase 2: 15 April 2020 – 3 May 2020 (19 days)
//     Phase 3: 4 May 2020 – 17 May 2020 (14 days)
//     Phase 4: 18 May 2020 – 31 May 2020 (14 days)
// Unlock:
//     Unlock 1.0: 1 June 2020 – 30 June 2020 (30 days)
//     Unlock 2.0: 1 July 2020 – 31 July 2020 (31 days)
//     Unlock 3.0: 1 August 2020 – 31 August 2020 (31 days)
//     Unlock 4.0: 1 September 2020 - 30 September 2020 (30 days)
//     Unlock 5.0: 1 October 2020 - 31 October 2020 (31 days)
//     Unlock 6.0: 1 November 2020 - 30 November 2020 (30 days)
//     Unlock 7.0: 1 December 2020 - 31 December 2020 (31 days)
//     Unlock 8.0: 1 January 2021 - 31 January 2021 (31 days)
//     Unlock 9.0: 1 February 2021 - 28 February 2021 (28 days)
//     Unlock 10.0: 1 March 2021 - 31 March 2021 (31 days)
//     Unlock 11.0: 1 April 2021 - 30 April 2021 (14 days)

  // Mar 2020 fit
  // if(t<15)  
  //  { sd = 0.7 + t*0.02/15.; }  // from 23rd Mar to 7th April
  // else if(t<36.) 
  //  { sd = 0.72 + (t-15.)*0.09/21.; }   // till 30th April 
  // else if(t<51.) 
  //  { sd = 0.81 - (t-36.)*0.06/15.; }  // till 15th May
  // else if(t<72.)
  //  { sd = 0.75 + (t-51.)*0.0125/21; }  // till 5 Jun  
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

   // Sep 2020 fit
      // if(t<30.)
      //  { sd = 0.77 + t*0.002/15; }
      // else if (t<153)
      //  { sd = 0.774 - (t-30)*0.003/15; } //till end of 2020
      // else
      //  { sd = 0.749 + (t-153)*delta/offsetsize; } 

  // sd = 0.7225; // no lockdown
  // Mar21: 15 days lockdown
  // if(t<35)  
  //   {  sd = 0.7225; }  // from 23rd Mar to 26th Apr,  // no lockdown
  // else if(t<50) 
  //  { sd = 0.7225 + (t-35.)*0.055/15.; }   // till   
  // else if(t<65.) 
  //  { sd = 0.7775 -(t-50.)*0.055/15.;  }
  // else
  //  { sd = 0.7225; }

  //Mar21: 21 days lockdown
  // if(t<35)  
  //   {  sd = 0.7225; }  // from 23rd Mar to 26th Apr,  // no lockdown
  // else if(t<56) 
  //  { sd = 0.7225 + (t-35.)*0.0775/21.; }   // till   
  // else if(t<71.) 
  //  { sd = 0.8 -(t-56.)*0.0775/15.;  } //unlock
  // else
  //  { sd = 0.7225; }

  //Mar21: 30 days lockdown (100%)
  // if(t<35)  
  //   {  sd = 0.7225; }  // from 23rd Mar to 26th Apr,  // no lockdown
  // else if(t<65) 
  //  { sd = 0.7225 + (t-35.)*0.116/30.; }   // till   
  // else if(t<80.) 
  //  { sd = 0.8385 -(t-65.)*0.116/15.;  } //unlock
  // else
  //  { sd = 0.7225; }
   

  //Mar21: 30 days lockdown  
  //updated on 20 May 2021, 14 day lockdown
  // if(t<35)  
  //   {  sd = 0.72; }  // from 23rd Mar to 26th Apr,  // no lockdown
  // else if(t<50) 
  //  { sd = 0.72 + (t-35.)*0.0348/15.; }   // till  
  // else if(t<60) 
  //  { sd = 0.7548 + (t-50.)*0.1352/10.; }   // till  
  // else if(t<70.) 
  //  { sd = 0.89 -(t-60.)*0.14/10.;  } //unlock, no locldown from 21st May
  // else
  //  { sd = 0.75; }

  // Mar21: 30 days lockdown  
  //updated on 20 May 2021, 14 day lockdown
  // if(t<35)  
  //   {  sd = 0.72; }  // from 23rd Mar to 26th Apr,  // no lockdown
  // else if(t<50) 
  //  { sd = 0.72 + (t-35.)*0.0348/15.; }   // till  
  // else if(t<60) 
  //  { sd = 0.7548 + (t-50.)*0.1352/10.; }   // till
  //  else if(t<75.) 
  //  { sd = 0.89;}  //extension of lockdown  for 15 days
  // else if(t<85.) 
  //  { sd = 0.89 -(t-75.)*0.14/10.;  } //unlock, no locldown from 21st May
  // else
  //  { sd = 0.75; }

  // Mar21: 30 days lockdown  
  //updated on 20 May 2021, 14 day lockdown
  // if(t<35)  
  //   {  sd = 0.72; }  // from 23rd Mar to 26th Apr,  // no lockdown
  // else if(t<50) 
  //  { sd = 0.72 + (t-35.)*0.0348/15.; }   // till  
  // else if(t<60) 
  //  { sd = 0.7548 + (t-50.)*0.1352/10.; }   // till
  //  else if(t<75.) 
  //  { sd = 0.89-(t-75.)*0.07/10.; }  //extension of lockdown for 15 days with 50%
  // else if(t<85.) 
  //  { sd = 0.82 -(t-75.)*0.07/10.;  } //unlock, no locldown from 21st May
  // else
  //  { sd = 0.75; }

  // Wave 3: Projections
  //updated on 04 Jun 2021 
  // if(t<15)  // no-age-wise distribution
  //   { sd = 0.6325; }  // from 23rd Mar to 26th Apr,  // no lockdown
  // else if(t<30) 
  //  { sd = 0.6325 + (t-15.)*0.005/15.; }   // till  
  // else
  //  { sd = 0.6375; }

  // three age groups' fit
  // if(t<42)  // from 23rd Mar to 28th Apr,  // no lockdown 
  //  { sd = 0.505; }  // day 35 lockdown + 7 days to start ralize
  // else if(t<72) // add lockdown effect in the next 15 days 
  //  { sd = 0.505 + (t-42.)*0.165/30.; }   // till  
  // else if(t<80)
  //  { sd = 0.67; } // 0.67 worked after peak
  // else if(t<90) // add lockdown effect in the next 15 days 
  //  { sd = 0.67 - (t-80.)*0.165/10.; }   // till  
  // else
  //  { sd = 0.57; }
    
  // //four age groups' fit    (wave II)
  // if(t<42)  // from 23rd Mar to 28th Apr,  // no lockdown 
  //  { sd = 0.455; }  // day 35 lockdown + 7 days to start ralize
  // else if(t<72) // add lockdown effect in the next 10 days 
  //  { sd = 0.455 + (t-42.)*0.215/30.; }   // till  
  // else if(t<80)
  //  { sd = 0.67; } // 0.67 worked after peak
  // else if(t<90) // add lockdown effect in the next 15 days 
  //  { sd = 0.67 - (t-80.)*0.215/10.; }   // till  
  // else
  //  { sd = 0.455; }

  //four age groups' fit Wave II, from 1 Jul 2020
  // if(t<30)  // from  
  //  { sd = 0.3; }  //  
  // if (t<120) //  
  //  { sd = 0.4 + (t-45.)*0.2/120.; }   // till  
  // else if(t<240) // 
  //  { sd = 0.6 - (t-120.)*0.2/30.; }   // till  
  // else 

  //fitted for CIR 1
  // if(t<30) // wave I - begin
  //  { sd = 0.525 ; } 
  // else if(t<40) 
  //  { sd = 0.525  + (t-30.)*0.045/10.; } 
  // else if(t<60)  
  //  { sd = 0.57; }    
  // else if(t<70)  
  //  { sd = 0.57  - (t-60.)*0.0225/10.; }     
  // else if(t<105)  
  //  { sd = 0.545; }       
  // else if(t<115) // wave I - ends
  //  { sd = 0.545  + (t-105.)*0.155/10.; }  
  // else if(t<125) { 
  //   sd = 0.7;}
  // else if(t<155) // Wave II - begin
  //  { sd = 0.7 - (t-125.)*0.175/30.; } 
  // else if(t<285) { 
  //   sd = 0.525; }     
  // else if(t<315) // Wave II - ends
  //  {sd = 0.525 + (t-285.)*0.0/30.; } 
  // else { 
  //   sd = 0.525;}

  // //fitted for CIR 1.75
  // if(t<30) // wave I - begin
  //  { sd = 0.515 ; } 
  // else if(t<35) 
  //  { sd = 0.515  + (t-30.)*0.04/5.; } 
  // else if(t<60)  
  //  { sd = 0.555; }   
  // else if(t<90)  
  //  { sd = 0.555  - (t-60.)*0.095/30.; }   
  // else if(t<115)  
  //  { sd = 0.46; }      //prev: 0.47 
  // else if(t<120) // wave I - ends
  //  { sd = 0.46  + (t-115.)*0.025/5.; }  
  // else if(t<310) // wave II - peak
  //  { sd = 0.485;} // prev: 0.48   
  // else if(t<330)  // wave I - begin
  //  { sd = 0.485 - (t-310.)*0.105/30.; } 
  // else { 
  //   sd = 0.38; } // prev: 0.4, (1L less) 

  // //fitted for CIR 1.25
  // if(t<30) // wave I - begin
  //  { sd = 0.52 ; } 
  // else if(t<35) 
  //  { sd = 0.52  + (t-30.)*0.045/5.; } 
  // else if(t<60)  
  //  { sd = 0.565; }   
  // else if(t<90)  
  //  { sd = 0.565 - (t-60.)*0.05/30.; }   
  // else if(t<100)  
  //  { sd = 0.515; }   
  // else if(t<105)  
  //  { sd = 0.515  + (t-100.)*0.095/5.; }  
  // else if(t<155) // wave II - peak
  //  { sd = 0.61;} // prev: .615 
  // else if(t<165) // tillmid Dec 2020
  //  { sd = 0.615 - (t-125.)*0.084/30.; } 
  // else { 
  //   sd = 0.521; } // prev: 0.512, (increase) :  0.53 (too high)

  // //fitted for CIR 1.25
  // if(t<30) // wave I - begin
  //  { sd = 0.525 ; } 
  // else if(t<35) 
  //  { sd = 0.525  + (t-30.)*0.03/5.; } 
  // else if(t<60)  
  //  { sd = 0.555; }   
  // else if(t<90)  
  //  { sd = 0.555 - (t-60.)*0.095/30.; }   
  // else if(t<130)  
  //  { sd = 0.46; }    // prev: 0.47 (slightly high, so reduced)
  // else if(t<135)  
  //  { sd = 0.46  + (t-130.)*0.048/5.; }  
  // // else if(t<155) // wave II - peak
  // //  { sd = 0.49;} // prev: .615 
  // // else if(t<165) // tillmid Dec 2020
  // //  { sd = 0.49 - (t-125.)*0.03/30.; } 
  // else { 
  //   sd = 0.508; } // prev: 0.49, (increase), 0.51, (decrease): 0.495 (increase): 0.50 (low, increase)
  //                 // 0.505(low, increase)
   
  
  if(TDatabase::ParamDB->KDP_nu==1)
   {
    // fitted for CIR 1.0, 9 Jul 2021, that is, 40 on an average
    if(t<25) // wave I - begin
     { sd = 0.71 ; } 
    else if(t<35) 
     { sd = 0.71  + (t-25.)*0.045/10.; } 
    else if(t<60)  
    { sd = 0.755; }  // prev: 0.745
    else if(t<70)  
    { sd = 0.755 - (t-60.)*0.031/10.; }   
    else if(t<100)  
     { sd = 0.721; }    // prev: 0.735 (reduced), 0.723 (high)
    else if(t<105)  
     { sd = 0.721  + (t-100.)*0.06/5; }  
     else if(t<150) 
    { sd = 0.781;}  // prev: .775 (low),  
    else if(t<220) // wave II - 
    { sd = 0.781 - (t-150.)*0.065/70.;} 
    else if(t<365) 
    { sd = 0.717; } //   0.721 (slightly high), 0.719 (low),  0.715 (too low), 0.719(high), 0.717(high)
    else
    { sd = 0.717 - TDatabase::ParamDB->KDP_lambda; }
   
   }
  else if(TDatabase::ParamDB->KDP_nu==1.25)
   {   
    //fitted for CIR 1.25, 9 Jul 2021, that is, 50 
    if(t<25) // wave I - begin
     { sd = 0.71 ; } 
    else if(t<35) 
     { sd = 0.71  + (t-25.)*0.03/10.; } 
    else if(t<60)  
    { sd = 0.74; }  // prev: 0.745
    else if(t<70)  
    { sd = 0.74 - (t-60.)*0.016/10.; }   
    else if(t<100)  
    { sd = 0.724; }    // prev: 0.73 (reduced), 0.725 (reduced), 0.722 (increase, 40k), 0.723 (increase, 40k)
    else if(t<105)  
    { sd = 0.724  + (t-100.)*0.026/5; }  
    else if(t<150) 
    { sd = 0.75;}  // prev: .775 (low), 0.74 (low) 
    else if(t<250) // wave II - 
    { sd = 0.75 - (t-150.)*0.043/100.;} 
    else if(t<365) 
    { sd = 0.707; }  // 0.71 (high, 80k),  0.70 (low, 1.1l),  0.705 (low, 25k),  0.708 (high)
    else { 
      sd = 0.707- TDatabase::ParamDB->KDP_lambda; } 
   }
  else
  {
   printf("TDatabase::ParamDB->KDP_nu = %f\n", TDatabase::ParamDB->KDP_nu);
#ifdef _MPI   
   MPI_Finalize();
#endif   
   exit(0); 
  }

  TDatabase::ParamDB->P13 =  (1. -  1./(1.+ exp(-(sd - sd_c)/bsd))); // f3 //fitted: 5.3* & P6 7.7
  // TDatabase::ParamDB->P13 = 4.;
  // 1. -  1./(1.+ exp(-( (0.3 + 0.125*sin(2.*22.*t/(30.*7.)) + (t-75.)*0.035/15.)   - 0.5)/0.1))
  //  OutPut("f1: " << TDatabase::ParamDB->P7 << " " );
  //  OutPut("f3: " << TDatabase::ParamDB->P13<<endl);
}

void GrowthAndB_Nuc_L0(int N_Inputs, double *Inn, double *Out) 
{
  double B_Nuc;
  double R_0  = TDatabase::ParamDB->P6; // Nucleation factor
  // double R_0  =3.35*1.01663385882;
  double f1  = TDatabase::ParamDB->P7; // interactive index
  // double bsigma  = TDatabase::ParamDB->P8;     
  // double sigma_c = TDatabase::ParamDB->P9;  
  double vsymp = TDatabase::ParamDB->P1;
  // double health  = TDatabase::ParamDB->P10;    
  // double bhealth = TDatabase::ParamDB->P11;  
  // double health_c = TDatabase::ParamDB->P12;  
  double sd  = TDatabase::ParamDB->P13; 
  double a4= 4.7, b4=0.28, c4_sqr=0.12*0.12;      
  // double t=TDatabase::TimeDB->CURRENTTIME ; 

  if(N_Inputs!=8)
  {
    printf("N_Inputs != 8,  Get_GrowthAndB_Nuc: %d \n",N_Inputs); 
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
  double VaccinePopRatio = Inn[6];
  double NormalisedVaccineRatio = Inn[7];  
  double VacFact = 0.05 + 0.95*VaccinePopRatio;

  double k5 = 2.3937;

  // cout << "VacFact: " <<  VacFact <<endl;
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

  // double f5, f4 = a4*exp(-(la -b4)* (la -b4)/c4_sqr); // used till 23 May 2021
  double f5, f4;
  
  //Jun 2021, https://ncdc.gov.in/dashboard.php
  // if(la<=10./125.) {
  //   f4 = 0.0336; }
  // else if(la<=20./125.) {
  //   f4 = 0.0841; }
  // else if(la<=30./125.){
  //   f4 = 0.218; }
  // else if(la<=40./125.){
  //   f4 = 0.2193; }    
  // else if(la<=50./125.){
  //   f4 = 0.1722; }    
  // else if(la<=60./125.){
  //   f4 = 0.1405; }    
  // else if(la<=70./125.){
  //   f4 = 0.0866; }    
  // else if(la<=80./125.){
  //   f4 = 0.0354; }    
  // else if(la<=90./125.){
  //   f4 = 0.0092; }    
  // else {
  //   f4 = 0.0012;}
  

 if (TDatabase::ParamDB->PBE_P10<0.)
  {
  if(la<=11./125.) {
    f4 = 2.56919365672*TDatabase::ParamDB->PBE_P5; 
    }
  else if(la<=17./125.) {
    f4 = 4.70419796961*TDatabase::ParamDB->PBE_P6; 
    }
  else if(la<=44./125.){
    f4 =  1.04537706409*TDatabase::ParamDB->PBE_P7;
     }
  else if(la<=59./125.){
    f4 = 1.88168560379*TDatabase::ParamDB->PBE_P8;
     }    
  else {
    f4 = 0.42765485306*TDatabase::ParamDB->PBE_P9;   
    }  
  f4 *=NormalisedVaccineRatio;      
  // cout << "NormalisedVaccineRatio: " << NormalisedVaccineRatio << endl;   
  }   
 else
 {
    f4 = NormalisedVaccineRatio/5.0; // N_AgeGroup = 5;
 }

  // 
  double lvsim2 = (lv-vsymp)*(lv-vsymp);
  if(lv<vsymp) {
    f5 = k5*exp( -lvsim2/(2.*(vsymp/3.)*(vsymp/3.)) ); }
   else  {
    f5 = k5*exp( -lvsim2/(VacFact*2*((1.-vsymp)/3.)*((1.-vsymp)/3.)) );  }
  
  double norm = TDatabase::ParamDB->FS_U;

  B_Nuc =SuscepPopRatio*R_0*f1*sd*f4*f5*IntL_val/norm; 

   if(B_Nuc<0)
   {
    cout << "B_Nuc: " << B_Nuc << " IntL_val "  <<IntL_val << " : " << sd << 
     " : " << f5 << endl;
    // exit(0);
   }  

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
  // G =  K_g*pow((la - arisk), p);
  G =  K_g* 0.114304; //mean value of pow((la - arisk), p);

  Out[0] = 10000000000000000; // growth
  Out[1] = 0;  //B_Nuc no nucleation
  Out[2] = 0; //L_Max Bound value,  if DIRICHLET 
   cout <<"Lv growth: " <<Out[0]<<endl;
   exit(0);
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

      //  if(lv<0.64) // mar 2020
       if(lv<0.64) //May 2021
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

        // Mar 21 
        // f1_Id = 0.007;
        // Mar 21 
       if(t<60)
         {         
          f1_Id = 0.02;
         }
       else if(t<70)
         { f1_Id = 0.02 - (t-60.)*0.012/10.;  }
       else if(t<270)      
         { f1_Id = 0.008; } // prev:  0.01;
       else if(t<290)
         { f1_Id = 0.008 + (t-270.)*0.003/20.;  }     
       else
         {
          f1_Id = 0.011; // 0.015 (high, reduce), 0.0125 (high, reduce)
         }


       } //  if(lv>=0.64)
        
      // Mar 2020 fit
      //  if(t<65)  //  till 26th May
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
      //   //  f1_R = 0.0475; // 17th June model
      //   f1_R = 0.0521 + (t-120.)*0.005/75; // 3 Aug model
      //   }  

        //  // Sep 2020 model
        //  if(t<30)
        //   {  f1_R =  0.04 + 1./(1+exp(-0.03*t))*0.75/15; }
        //  else
        //   { f1_R =  0.076; } 
        
      //   //Mar/May 2021
      //  if(t<30)
      //    { 
      //     f1_R =  0.04 + t*0.0075/30.;
      //    }
      //   else if(t<60) 
      //    { 
      //     f1_R = 0.0475 + (t-30.)*0.0235/30;  
      //    }
      //  else
      //   {
      //    f1_R = 0.071;
      //   }   

       //   if(t<40) // policy change on 27th Apr 2021 : Wave II fit
        //    {     
        //     f1_R = 0.05;
        //    }         
        //   else if(t<50) 
        //   { 
        //    f1_R = 0.05 + (t-40.)*0.02/10;  
        //   }
        //  else
        //   { 
        //    f1_R = 0.07;
        //   }         

           // initially recovery has to be small, otherwise, both acive and recovery will be very low
 
         if(t<240)  
           {     
            f1_R = 0.075;
           }         
          else if(t<250) 
          { 
           f1_R = 0.075 - (t-240.)*0.0025/10;  
          }
         else  
          { 
           f1_R = 0.0725;
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
                      DoubleFunctND **GrowthAndB_Nuc, double **XPos)
  {
   // L_d  
   L_StartEnd[0] = TDatabase::ParamDB->REACTOR_P12;
   L_StartEnd[1] = TDatabase::ParamDB->REACTOR_P13;
   L_NVertices[0] = (int)TDatabase::ParamDB->REACTOR_P11;
   GrowthAndB_Nuc[0] = GrowthAndB_Nuc_L0;
   BoundConLminLMax[0] = BoundCondition_LminLMax_L0;
   XPos[0] = nullptr;

   // L_v  
   L_StartEnd[2] = 0.;
   L_StartEnd[3] = 1.;
   L_NVertices[1] = (int)TDatabase::ParamDB->REACTOR_P15;
   GrowthAndB_Nuc[1] = GrowthAndB_Nuc_L1;
   BoundConLminLMax[1] = BoundCondition_LminLMax_L1;
   XPos[1] = nullptr;

   // L_a  
   L_StartEnd[4] = 0.;
   L_StartEnd[5] = 1.0; // la/125
   L_NVertices[2] = (int)TDatabase::ParamDB->REACTOR_P16;
   GrowthAndB_Nuc[2] = GrowthAndB_Nuc_L2;
   BoundConLminLMax[2] = BoundCondition_LminLMax_L2;

  //  double X_LA[] = {0., 5.5/125., 11./125., 14/125., 17./125., 30.5/125., 44./125., 
  //                  51.5/125., 59./125., 92./125, 125./125}; 
   double X_LA[] = {0., 11./125., 17./125., 44./125., 59./125., 125./125}; 
   L_NVertices[2] = 6;

  for(int i=0; i< L_NVertices[2]; i++)
    XPos[2][i] = X_LA[i];

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
