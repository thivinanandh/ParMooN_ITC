// =======================================================================
// Purpose:     header for Covid model with new kernels of ParMooN
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 04.04.2020
// =======================================================================

#include <ADISystem.h>
#include <ADISystem1D.h>
#include <LineAffin.h>
#include <NodalFunctional1D.h>
#include <QuadIsoparametric.h>
#include <QuadBilinear.h>
#include <TriaIsoparametric.h>
#include <MooNMD_Io.h>

#include <MacroCell.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
// #include <TimeConvDiff2D.h>
#include <MacroCell.h>
#include <iostream>
using namespace std;

#include "StatePopulationData.h"

void WriteInitAntiBody(double Value, char *Name)
{
 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
//  os << "PopulationData/Covid_Antibody.data" << ends;
 os << "PopulationData/Covid_"<<Name<<".data" << ends;
 std::ofstream dat(os.str().c_str());
 dat.setf(std::ios::fixed);
 
 if (!dat)
  {
   cerr << "cannot open file for output" << endl;
   exit(0);
  }

  dat << "Covid-19 Antibody Estimated by PSD Model @ CDS/IISc" << endl;
  dat << "Day, AntiBody" <<endl;
  // cout << "N_AllCell : " << N <<endl;
  dat << setprecision(0) << TDatabase::TimeDB->CURRENTTIME;
  dat <<" "<< Value <<endl;

  dat.close();
//  OutPut( "Covid "<<Name<<" data wrote into file " <<endl);
}

void WriteAntiBody(double Value, char *Name)
{
 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
//  os << "PopulationData/Covid_Antibody.data" << ends;
 os << "PopulationData/Covid_"<<Name<<".data" << ends;

 std::ofstream dat(os.str().c_str(), std::ios_base::app);
 dat.setf(std::ios::fixed);
  
 dat << setprecision(0) << TDatabase::TimeDB->CURRENTTIME;
 dat << " " << Value <<endl;

 dat.close();
} // WriteData

void WriteInitData(int N, double *Values, char *Name)
{
 double sum=0;
 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
 os << "PopulationData/Covid"<<Name<<".data" << ends;
  
 std::ofstream dat(os.str().c_str());
 dat.setf(std::ios::fixed);
 
 if (!dat)
  {
   cerr << "cannot open file for output" << endl;
   exit(0);
  }

  dat << "Covid-19 "<<Name<<" Estimated by PSD Model @ CDS/IISc" << endl;
  dat << "Day" ;
  // cout << "N_AllCell : " << N <<endl;

 for(int i=0; i<N; ++i)
   dat<< " " << enum_to_string(i);
 dat << " Total" <<endl;

 dat << setprecision(1) << TDatabase::TimeDB->CURRENTTIME;
 for(int i=0; i<N; ++i)
  { 
   dat << " " << setprecision(0) << Values[i];
   sum +=Values[i];
  }
 dat << " " << sum <<endl;

 dat.close();
 cout << endl;
 OutPut( "Covid "<<Name<<" data wrote into file " <<endl);
} // WriteData

void WriteData(int N, double *Values, char *Name)
{
 double sum;
 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
 os << "PopulationData/Covid"<<Name<<".data" << ends;
 
 std::ofstream dat(os.str().c_str(), std::ios_base::app);
 dat.setf(std::ios::fixed);
  
 sum=0.;
 dat << setprecision(1) << TDatabase::TimeDB->CURRENTTIME;
 for(int i=0; i<N; ++i)
  { 
   dat << " " << setprecision(0) << Values[i];
  //  cout<<"Sum " << sum <<endl;
   sum +=Values[i];
  }
 dat << " " << sum <<endl;

 dat.close();
//  cout << endl;
//  OutPut( "Covid Population data wrote into file " << N_IData <<endl);
 
} // WriteData

void WriteInitDistData(int *N_LnodalPos, double **LnodalPos, int N_LDistXSum, double *LDistXSum)
{
 int i; 
 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
 os << "PopulationData/CovidDist.data" << ends;
  
 std::ofstream dat(os.str().c_str());
 dat.setf(std::ios::fixed);
 
 if (!dat)
  {
   cerr << "cannot open file for output" << endl;
   exit(0);
  }

  dat << "Covid-19 Distribution Estimated by PSD Model @ CDS/IISc. First three rows contain l_a, l_v, l_d positions, subsequent row contain the day and Active data " << endl;

 for(i=0; i<N_LnodalPos[2]-1; ++i)
   dat<< " " << setprecision(2) << LnodalPos[2][i]<< ", " ;

 dat<< " " << setprecision(2) << LnodalPos[2][N_LnodalPos[2]-1];
 dat <<endl;

 for(i=0; i<N_LnodalPos[1]-1; ++i)
   dat<< " " << setprecision(2) <<  LnodalPos[1][i]<< ", ";
 
 dat<< " " << setprecision(2) <<  LnodalPos[1][N_LnodalPos[1]-1];   
 dat <<endl;

 for(i=0; i<N_LnodalPos[0]-2; i=i+2)
   dat<< " " << setprecision(2) << (LnodalPos[0][i]+LnodalPos[0][i+1])/2.<< ", "; // write sol only for mid point, ANSATZ_ORDER_INTL: -11 

 dat<< " " << setprecision(2) << (LnodalPos[0][N_LnodalPos[0]-2]+LnodalPos[0][N_LnodalPos[0]-1])/2.; // write sol only for mid point, ANSATZ_ORDER_INTL: -11 
 dat <<endl;


 dat << setprecision(1) << TDatabase::TimeDB->CURRENTTIME << ", ";
 for(int i=0; i<N_LDistXSum-1; ++i)
  dat << " " << setprecision(4) << LDistXSum[i]<< ", ";

 dat << " " << setprecision(4) << LDistXSum[N_LDistXSum-1];  
  dat <<endl;

 dat.close(); 
 OutPut( "Covid CovidDist data wrote into file " <<endl);
}


void WriteDistData(int N_LDistXSum, double *LDistXSum)
{
 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
 os << "PopulationData/CovidDist.data" << ends;

 std::ofstream dat(os.str().c_str(), std::ios_base::app);
 dat.setf(std::ios::fixed);

 dat << setprecision(1) << TDatabase::TimeDB->CURRENTTIME<< ", ";
 for(int i=0; i<N_LDistXSum-1; ++i)
  dat << " " << setprecision(4) << LDistXSum[i]<< ", ";

 dat << " " << setprecision(4) << LDistXSum[N_LDistXSum-1];  
 dat <<endl;

 dat.close(); 
//  OutPut( "Covid CovidDist data wrote into file " <<endl);
}



void WriteInitRecovDistData(int *N_LnodalPos, double **LnodalPos, int N_LDistXSum, double *LRecoDistXSum)
{
 int i; 
 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
 os << "PopulationData/CovidRecovDist.data" << ends;
  
 std::ofstream dat(os.str().c_str());
 dat.setf(std::ios::fixed);
 
 if (!dat)
  {
   cerr << "cannot open file for output" << endl;
   exit(0);
  }

  dat << "Covid-19 Distribution Estimated by PSD Model @ CDS/IISc. First three rows contain l_a, l_v, l_d positions, subsequent row contain the day and recoverd data " << endl;

 for(i=0; i<N_LnodalPos[2]-1; ++i)
   dat<< " " << setprecision(2) << LnodalPos[2][i]<< ", " ;

 dat<< " " << setprecision(2) << LnodalPos[2][N_LnodalPos[2]-1];
 dat <<endl;

 for(i=0; i<N_LnodalPos[1]-1; ++i)
   dat<< " " << setprecision(2) <<  LnodalPos[1][i]<< ", ";
 
 dat<< " " << setprecision(2) <<  LnodalPos[1][N_LnodalPos[1]-1];   
 dat <<endl;

 for(i=0; i<N_LnodalPos[0]-2; i=i+2)
   dat<< " " << setprecision(2) << (LnodalPos[0][i]+LnodalPos[0][i+1])/2.<< ", "; // write sol only for mid point, ANSATZ_ORDER_INTL: -11 

 dat<< " " << setprecision(2) << (LnodalPos[0][N_LnodalPos[0]-2]+LnodalPos[0][N_LnodalPos[0]-1])/2.; // write sol only for mid point, ANSATZ_ORDER_INTL: -11 
 dat <<endl;


 dat << setprecision(1) << TDatabase::TimeDB->CURRENTTIME << ", ";
 for(int i=0; i<N_LDistXSum-1; ++i)
  dat << " " << setprecision(4) << LRecoDistXSum[i]<< ", ";

 dat << " " << setprecision(4) << LRecoDistXSum[N_LDistXSum-1];  
  dat <<endl;

 dat.close(); 
 OutPut( "Covid CovidRecovDist data wrote into file " <<endl);
}

void WriteRecovDistData(int N_LDistXSum, double *LRecoDistXSum)
{
 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
 os << "PopulationData/CovidRecovDist.data" << ends;

 std::ofstream dat(os.str().c_str(), std::ios_base::app);
 dat.setf(std::ios::fixed);

 dat << setprecision(1) << TDatabase::TimeDB->CURRENTTIME<< ", ";
 for(int i=0; i<N_LDistXSum-1; ++i)
  dat << " " << setprecision(4) << LRecoDistXSum[i]<< ", ";

 dat << " " << setprecision(4) << LRecoDistXSum[N_LDistXSum-1];  
 dat <<endl;

 dat.close(); 
//  OutPut( "Covid CovidDist data wrote into file " <<endl);
}

void ReadData(int N, char *DataFile, double *PopulationData)
{
 char line[100], aux_char[251];
 double val, val1=0., total=0.;
 string temp;
 int n=0;

 std::ifstream dat(DataFile);

 if (!dat)
  {
   cerr << "cannot open file: " << DataFile << endl;
   exit(0);
  }

 //  memset(PopulationData, 0, 32*SizeOfInt);
 for(int i=0; i<N; i++)
  {
   dat >> line;
    
   temp= enum_to_string(i);   
   strcpy(aux_char, temp.c_str()); 
 
   //code to generate switch code
  //  cout<< "case " << i <<":"<<endl;
  //     cout << "    return \""<<line <<"\";"<< endl;
   // cout<< line <<", ";
 
   if(!strcmp(line, aux_char))
    {
     dat >> val;
     PopulationData[i]= val;  
     total +=val;
    // if(TDatabase::ParamDB->Par_P0==1)   
    //  {    
    //   val1 +=val;          
    //   cout <<val1 <<" Population: " <<  aux_char <<  " " << PopulationData[i]<< endl;
    //  }
    }
    // read until end of line
    dat.getline (line, 99);
  } // for i
  // cout<< "Total " << total<<endl;
 dat.close(); 
}

void ReadPopulationData(int N, char *DataFile, int N_Days, double **AntibodyPopulation, double fact)
{
 char line[1080], aux_char[20];
 double val, val1;
 string temp;
 int j, n=0;
 std::string str; 

 std::ifstream dat(DataFile);

 if (!dat)
  {
   cerr << "cannot open file: " << DataFile << endl;
   exit(0);
  }

 //  memset(PopulationData, 0, 32*SizeOfInt);
 for(int i=0; i<N; i++)
  {
   dat >> line;
    
   temp= enum_to_string(i);   
   strcpy(aux_char, temp.c_str()); 
 
   //code to generate switch code
  //   if(TDatabase::ParamDB->Par_P0==1) 
  //   {
  // //  cout<< "temp " <<temp <<" : " << line <<endl;
  //     // cout << "    return \""<<line <<"\";"<< endl;
  // //  cout<< line <<", ";
  //   }
   if(!strcmp(line, aux_char))
    {
      for(j=0; j<N_Days; j++)
       {
        dat >> val;
        AntibodyPopulation[i][j]= fact*val;  

        // if(TDatabase::ParamDB->Par_P0==1)           
        // cout << i << " " << j << " " << val << " Population: " <<  aux_char <<  " " << AntibodyPopulation[i][j]<< endl;
              
       } // for(j=0; j<250; j++)
    }
// MPI_Finalize();
// exit(0);
    // read until end of line
    dat.getline (line, 250);
  } // for i
  // cout<<endl;
 dat.close(); 
}

void Sort(double *XArray, double *YArray, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  double Mid, Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=XArray[m];

      do
      {
        while( XArray[i] > Mid) i++;

        while(XArray[j] < Mid) j--;

        if (i<=j)
        {
          Temp=XArray[i];
          XArray[i]=XArray[j];
          XArray[j]=Temp;

          Temp=YArray[i];
          YArray[i]=YArray[j];
          YArray[j]=Temp;

          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete [] rr;

}


int GetIndex(double *x, double *y, int Length, double x0, double y0)
{
  int l=0, r=Length, m=(r+l)/2;
  double MidX, MidY;
  bool update;

  MidX=x[m];
  MidY=y[m];

  //   cout<<endl;
  //   cout << " x " << x0<< " y " << y0<<endl;
  //   cout << " MidX " << MidX<< " MidY " << MidY<<endl;
  //  if((fabs(MidX-x0)<1e-8 &&  fabs(MidY-y0)<1e-8))
  //     cout << "TRUE " <<endl;

  while(!( fabs(MidX-x0)<1.e-4 && fabs(MidY-y0)<1.e-4 ) )
  {

    if(fabs(MidX-x0) <1.e-4 )
    {
      if(MidY> y0)
        { l=m;}
        else
          { r=m; }
    }
    else if(MidX > x0)
      { l=m; }
      else
        { r=m; }

        m=(r+l)/2;
    MidX=x[m];
    MidY=y[m];
    //cout << m << " MidX " << MidX<< " MidY " << MidY<< " x0  " <<  x0 << "  y0 " <<  y0 <<endl;

  }
  //  cout << " outer m " << m<< endl;
  //   cout << r << " MidX " << MidX<< " MidY " << MidY<< " y[m-1] " <<  y[m-1]<< "  y[m+1] " <<  y[m+1]<<endl;
  return m;
}

int GetIndex(double *x, int Length, double x0, double y0)
{
  int l=0, r=Length, m=(r+l)/2;
  double MidX, MidY;
  bool update;

  MidX=x[2*m];
  MidY=x[2*m+1];

  //   cout<<endl;
  //   cout << " x " << x0<< " y " << y0<<endl;
  //   cout << " MidX " << MidX<< " MidY " << MidY<<endl;
  //  if((fabs(MidX-x0)<1e-8 &&  fabs(MidY-y0)<1e-8))
  //     cout << "TRUE " <<endl;

  while(!( fabs(MidX-x0)<1.e-4 && fabs(MidY-y0)<1.e-4 ) )
  {

    if(fabs(MidX-x0) <1.e-4 )
    {
      if(MidY> y0)
        { l=m;}
        else
          { r=m; }
    }
    else if(MidX > x0)
      { l=m; }
      else
        { r=m; }

        m=(r+l)/2;
    MidX=x[2*m];
    MidY=x[2*m+1];
    //cout << m << " MidX " << MidX<< " MidY " << MidY<< " x0  " <<  x0 << "  y0 " <<  y0 <<endl;

  }
  //  cout << " outer m " << m<< endl;
  //   cout << r << " MidX " << MidX<< " MidY " << MidY<< " y[m-1] " <<  y[m-1]<< "  y[m+1] " <<  y[m+1]<<endl;
  return m;
}


void GetNodalPts(int &N_Xpos, TFESpace2D *FESpace2D, int &MaxN_PtsForNodal, double *&IntlX)
{
  int i,j,k,l, N_AllLocalPoints, N_V;
  int N_Cells, PolynomialDegree;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, *RootPtIndex;
  int *DOF, N_Edges, ApproxOrder, N_XCommonPts, N_RootPts, disp;

  double *xi, *eta, x0, y0;
  double *x_loc, *y_loc, *x_loc_origOrder, *y_loc_origOrder;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D], lastx, lasty;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];

  bool IsIsoparametric, update;

  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TFE2D *FE_Obj;
  TNodalFunctional2D *nf;
  TBaseFunct2D *bf;
  RefTrans2D RefTrans, *RefTransArray;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  BF2DRefElements RefElement;
  TRefTrans2D *rt;
  RefTrans2D F_K;
  QuadFormula2D QuadFormula;
  TVertex *Vertices;

  Coll = FESpace2D->GetCollection();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_Cells = Coll->GetN_Cells();

  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

  //first find the number of AllLocal points
  N_AllLocalPoints = 0;
  MaxN_PtsForNodal = -1;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_AllLocalPoints +=N_Points;
    if(MaxN_PtsForNodal<N_Points) MaxN_PtsForNodal=N_Points;
  }

  x_loc = new double[N_AllLocalPoints];
  y_loc = new double[N_AllLocalPoints];
  x_loc_origOrder = new double[N_AllLocalPoints];
  y_loc_origOrder = new double[N_AllLocalPoints];

  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);

    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ApproxOrder = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    RefElement = Element->GetBaseFunct2D()->GetRefElement();

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        N_Edges = 4;
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        N_Edges = 3;
        break;
    }

    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint ||
          jointtype == IsoBoundEdge)
          IsIsoparametric = TRUE;
      }
    }                                             // endif

    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          break;
      }
    }

    switch(RefTrans)
    {
      case QuadAffin:
        rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
        ((TQuadAffin *)rt)->SetCell(cell);
        F_K = QuadAffin;
        break;
      case QuadBilinear:
        rt = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)rt)->SetCell(cell);
        F_K = QuadBilinear;
        break;
      case QuadIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(QuadIsoparametric);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        F_K = QuadIsoparametric;
        break;
      case TriaAffin:
        rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)rt)->SetCell(cell);
        F_K = TriaAffin;
        break;
      case TriaIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(TriaIsoparametric);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        F_K = TriaIsoparametric;
        break;
      default:
        cout << "unknown reftrans id: " << RefTrans << endl;
    }

    TFEDatabase2D::GetOrigFromRef(F_K, N_Points, xi, eta,  X, Y, AbsDetjk);

    for(j=0; j<N_Points; j++)
    {
      x_loc[N_AllLocalPoints] = X[j];
      y_loc[N_AllLocalPoints] = Y[j];
      N_AllLocalPoints++;
    }
  }                                               // for(i=0; i<N_Cells; i++)

  memcpy(x_loc_origOrder, x_loc,  N_AllLocalPoints*SizeOfDouble);
  memcpy(y_loc_origOrder, y_loc,  N_AllLocalPoints*SizeOfDouble);


  // reove common/shared nodal points 
  Sort(x_loc, y_loc, N_AllLocalPoints);

  lastx  = x_loc[0];
  N_RootPts = 0;
  N_XCommonPts = 0;
  disp = 0;

  for(i=0; i<N_AllLocalPoints; i++)
  {
    if( fabs(x_loc[i]-lastx)>1e-5 )
    {
      xi = x_loc+disp;
      eta = y_loc+disp;
      Sort(eta, xi, N_XCommonPts);

      lasty= eta[0];
      N_RootPts++;
      for(j=0; j<N_XCommonPts; j++)
      {
        if( fabs(eta[j]-lasty)>1e-4)
        {
          N_RootPts++;
          lasty = eta[j];
        }
      }                                           //  for(j=0; j<N_XCommonPts; j++)

      disp +=N_XCommonPts;
      N_XCommonPts = 1;
      lastx = x_loc[i];
    }
    else
      { N_XCommonPts++; }
    }                                             // for(i=0; i<N_AllLocalPoints; i++)
    // sort the final part
    xi = x_loc+disp;
  eta = y_loc+disp;
  Sort(eta, xi, N_XCommonPts);
  lasty= eta[0];
  N_RootPts++;
  for(j=0; j<N_XCommonPts; j++)
  {
    if( fabs(eta[j]-lasty)>1e-4)
    {
      N_RootPts++;
      lasty = eta[j];
    }
  }                                               //  for(j=0; j<N_XCommonPts; j++)

  // cout<< "N_RootPts " << N_RootPts <<endl;
  // exit(0);
  IntlX = new double[2*N_RootPts];
  // IntlY = new double[N_RootPts];

  lastx  = x_loc[0];
  N_RootPts = 0;
  N_XCommonPts = 0;
  disp = 0;
  for(i=0; i<N_AllLocalPoints; i++)
  {
    if( fabs(x_loc[i]-lastx)>1e-4 )
    {
      xi = x_loc+disp;
      eta = y_loc+disp;
      lasty= eta[0];
      IntlX[2*N_RootPts] = xi[0];
      IntlX[2*N_RootPts+1] = eta[0];
      N_RootPts++;

      for(j=0; j<N_XCommonPts; j++)
      {
        if( fabs(eta[j]-lasty)>1e-4)
        {
          IntlX[2*N_RootPts] = xi[j];
          IntlX[2*N_RootPts+1] = eta[j];
          N_RootPts++;
          lasty = eta[j];
        }
      }                                           //  for(j=0; j<N_XCommonPts; j++)
      disp +=N_XCommonPts;
      N_XCommonPts = 1;
      lastx = x_loc[i];
    }
    else
      { N_XCommonPts++; }
    }                                             //for(i=0; i<N_AllLocalPoints; i++)

    xi = x_loc+disp;
  eta = y_loc+disp;
  lasty= eta[0];
  IntlX[2*N_RootPts] = xi[0];
  IntlX[2*N_RootPts+1] = eta[0];
  N_RootPts++;
  for(j=1; j<N_XCommonPts; j++)
  {
    if( fabs(eta[j]-lasty)>1e-4)
    {
      IntlX[2*N_RootPts] = xi[j];
      IntlX[2*N_RootPts+1] = eta[j];
      N_RootPts++;
      lasty = eta[j];
    }
  }                                               //  for(j=0; j<N_XCommonPts; j++)

  //   for(i=0; i<N_RootPts; i++)
  //    cout<< i << " IntlX " << IntlX[2*i]<< " IntlY " << IntlY[2*i+1]  <<endl;

  delete [] x_loc;
  delete [] y_loc;

  RootPtIndex = new int[N_AllLocalPoints];

  // find the index for the local points in the root points
  for(i=0; i<N_AllLocalPoints; i++)
  {
    x0 = x_loc_origOrder[i];
    y0 = y_loc_origOrder[i];

    RootPtIndex[i]=GetIndex(IntlX, N_RootPts,  x0, y0);
  }

  delete [] x_loc_origOrder;
  delete [] y_loc_origOrder;

  FESpace2D->SetIntlPtIndexOfPts(RootPtIndex);

  N_Xpos = N_RootPts;

  //    for(i=0; i<N_Xpos; i++)
  //        cout <<i<< " IntlX[i] " <<IntlX[2*i]<< " IntlY[2*i+1] " <<IntlY[i] <<endl;
  
  // cout<< "N_RootPts " << N_RootPts <<endl;
  // exit(0);
}



void GetQuadPts(int &N_Xpos, TFESpace2D *FESpace2D, double *&IntlX)
{
  int i,j,k,l, m;
  int N_Cells, PolynomialDegree, N_Points;

  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];
  double *xi, *eta, *weights;

  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  TRefTrans2D *F_K;

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  m = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)F_K)->SetCell(cell);
        break;
    }                                             // endswitch

    //cout << "QuadFormula: " << QuadFormula << endl;
    qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
    qf2->GetFormulaData(N_Points, weights, xi, eta);

    if(i==0)
    {
      N_Xpos = N_Points*N_Cells;
      //     cout<< "N_Xpos " << N_Xpos << endl;
      IntlX = new double[2*N_Xpos];
      // IntlY = new double[N_Xpos];
    }

    switch(RefElement)
    {
      case BFUnitSquare:
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;

      case BFUnitTriangle:
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;
    }                                             // endswitch

    for(j=0; j<N_Points; j++)
    {
      IntlX[2*m] = X[j];
      IntlX[2*m+1] = Y[j];
      m++;
    }
  }                                               // for(i=0; i<N_Cells; i++)

  //  cout<< "N_Xpos " << m <<endl;
  //   for(i=0; i<N_Xpos; i++)
  //    cout<< i << " IntlX " << IntlX[2*i]<< " IntlY " << IntlX[2*i+1]  <<endl;
  //  exit(0);

}

void WeibulDistribution(int N, double *WeibulDist)
{
  double k = 5.67, lamb;  
  // double k = 5.67, lamb = 150.; // Siva's paper lamb = 120

  if(TDatabase::TimeDB->CURRENTTIME<365.)
  { lamb = 180.; }
  else
  { lamb = TDatabase::ParamDB->KDP_g_2; }

  // Weibul Distribution
  WeibulDist[0] = 0.;
  for(int i=1; i<N; i++)
  {
   WeibulDist[i] += WeibulDist[i-1] 
                    + (k/lamb)*pow(double(i)/lamb, (k-1.))*exp(-pow(double(i)/lamb, k));                    
  }

  //cummulative
  for(int i=0; i<N; i++)
  {
   WeibulDist[i] =1. -  WeibulDist[i];
  //  cout << i << " WeibulDist[i] " << WeibulDist[i] <<endl; 
  }
}

void DistributeWaningAntibody(int N_Dist, int Today, int N, double *Val, double *WeibulDist, 
                                double **AntibodyPopulation, double ratio, 
                                double *State_CIR, double CIR, int NewVariant, double *StatePopulation)
  {
   int i, j, l;
   double val_dist=1., value;
   double k = 2.105;
   double lamb = 40.;   
   double CIRFactor = TDatabase::ParamDB->KDP_nu;
   double ABWPop = TDatabase::ParamDB->KDP_g_1;

  if(Today<365)
  { ABWPop = 1.; }
  else
  { ABWPop = TDatabase::ParamDB->KDP_g_1; }


   for(j=Today; j<N; j++)
      {
       l=j-Today;
       if(NewVariant>0)
        val_dist = 1. - 45.*(k/lamb)*pow(double(NewVariant+l)/lamb, (k-1.))*exp(-pow(double(NewVariant+l)/lamb, k));

       value = val_dist*WeibulDist[l];

       for(i=0; i<N_Dist; i++)
        {
         if(CIR<0) {
          CIR = State_CIR[i]*CIRFactor; }
          
        //  AntibodyPopulation[i][j] += value*Val[i]*CIR;
        AntibodyPopulation[i][j] += ((1.-ABWPop) +  ABWPop*value)*Val[i]*CIR;

        if(AntibodyPopulation[i][j] > ratio*StatePopulation[i])
         AntibodyPopulation[i][j] = ratio*StatePopulation[i]; 

      } // i
  } // j

// cout <<"maxTPR: " << maxTPR <<endl;
// MPI_Finalize();
// exit(0);
  }

 //waning of existing antibody
 void InitNewVirusVariant(int today, int N_AgeGroup, int N_Xpos, int N_EffAntibdyDays, 
                        double ***AntibodyPopulation)
 {
  // Weibul distribution
  double val=0., k = 2.67;
  double lamb = 15.;
  int i, j, l;
  // int end = today + N_EffAntibdyDays;
  int end = int(TDatabase::TimeDB->ENDTIME);

  double ImmuneEscapeVal;
  
  if(today<365)
  { ImmuneEscapeVal = 1.; }
  else
  { ImmuneEscapeVal = TDatabase::ParamDB->KDP_w_sat_1; }

  for(i = today; i <end; i++)
   {
    val += (k/lamb)*pow(double(i-today)/lamb, (k-1.))*exp(-pow(double(i-today)/lamb, k));

    for(j=0; j<N_AgeGroup; j++)
     for(l=0; l<N_Xpos; l++)
      {
      //  AntibodyPopulation[j][l][i]  *=(1. - val);
       AntibodyPopulation[j][l][i]  = ((1.-ImmuneEscapeVal) + ImmuneEscapeVal*(1. - val))*AntibodyPopulation[j][l][i] ;      
      //  if(TDatabase::ParamDB->SDFEM_TYPE==0 && j==0 && l==0)
      //   cout<<i << " AntibodyPopulation[j][l][i]: " << AntibodyPopulation[j][l][i] << " " <<  1. - val 
      //     << " " << AntibodyPopulation[j][l][i]*(1. - val) << endl;

      }

    // cout<<i << " NewVarFact: " << val << endl;
   }
 }
