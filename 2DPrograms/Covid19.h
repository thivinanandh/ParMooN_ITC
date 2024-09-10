

#include <ADISystem.h>
#include <ADISystem1D.h>
#include <LineAffin.h>
#include <NodalFunctional1D.h>
#include <QuadIsoparametric.h>
#include <QuadBilinear.h>
#include <TriaIsoparametric.h>
#include <MooNMD_Io.h>

 

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


void Generate1DMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, *X;
  double hmin, hmax;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  N_V = N_Cells+1;
  X = new double[N_V];

  h = (End-Start)/(double)N_Cells;

  X[0] = Start;

  for(i=1; i<N_V; i++)
   X[i] = X[i-1] + (double)h;

  X[N_V-1] = End;

    hmin = 1.e8; hmax = -1.e8; 
    for(i=0; i<N_V-1; i++)
     {
      len = sqrt ((X[i+1] - X[i])*(X[i+1] - X[i]));
      if(len< hmin) hmin = len;
      if(len> hmax) hmax = len;        
     }
     OutPut("L h_min : " << hmin << " L h_max : " << hmax << endl);
    
 //   for(i=0; i<N_V; i++)
 //    cout<< i << " X[i] " << X[i] <<endl;
 // 
 // exit(0);

  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_V]; 

  y=0.;

  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
    Vetrex[i] = new TVertex(X[i], y);
   }

  Vetrex[N_Cells] = new TVertex(X[N_V-1], y);

  CellTree = new TBaseCell*[N_Cells];

   for (i=0;i<N_Cells;i++)
   {
 //     Vetrex[ i ]->GetCoords(x, y);
 //     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);
    ((TMacroCell *) CellTree[i])->SetSubGridID(0);
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


   for(i=1;i<N_Cells;i++)
    {
     Joint = new TJointEqN(CellTree[i-1], CellTree[i]);

     CellTree[i-1]->SetJoint(1, Joint);
     CellTree[i]->SetJoint(0, Joint);
   } // for(i=0;i<N_Cells;i++)

   // end joint(vertex)
   Joint = new TJointEqN(CellTree[N_Cells-1]);
   CellTree[N_Cells-1]->SetJoint(1, Joint);

  delete []  Lines;
}


void GetInternalNodalPts(TFESpace1D * FeSpace_L1, int &N_L1nodal, double *&IntlPosL)
{
  int i, j, k, l, m, r, N_Cells, N_RootPts, *RootPtIndex;
  int N_DOFs, N_LocalDOFs, N_Points;
  int N_AllLocalPoints;

  double L, L0, *xi, *eta, *L_loc, *L_loc_origOrder;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;

  Coll = FeSpace_L1->GetCollection();
  N_Cells = Coll->GetN_Cells();

  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FeSpace_L1->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    N_AllLocalPoints +=N_Points;
  }

  L_loc = new double [N_AllLocalPoints];
  L_loc_origOrder = new double [N_AllLocalPoints];
  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FeSpace_L1->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);

    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, X, Y, AbsDetjk);

    for(j=0; j<N_Points; j++)
    {
      L_loc[N_AllLocalPoints] = X[j];
      N_AllLocalPoints++;
    }
  } // for(i=0; i<N_Cells; i++)

  memcpy(L_loc_origOrder, L_loc,  N_AllLocalPoints*SizeOfDouble);

  for(i=0; i<N_AllLocalPoints-1; i++)
   for(j=i; j<N_AllLocalPoints; j++)
    if(L_loc[i]> L_loc[j])
     {
      L= L_loc[i];
      L_loc[i]= L_loc[j];
      L_loc[j]= L;
     }

  L  = L_loc[0];
  N_RootPts = 1;

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      N_RootPts++;
      L = L_loc[i];
     }
   }

  IntlPosL= new double[N_AllLocalPoints];
  IntlPosL[0] = L_loc[0];
  N_RootPts = 1;
  L  = L_loc[0];

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      IntlPosL[N_RootPts] = L_loc[i];
      N_RootPts++;
      L = L_loc[i];
     }
   }

  delete [] L_loc;

  RootPtIndex = new int[N_AllLocalPoints];

  // find the index for the local points in the root points
  for(i=0; i<N_AllLocalPoints; i++)
   {
    L = L_loc_origOrder[i];
    l=0;
    r=N_RootPts;

    m = N_RootPts/2;
    L0 = IntlPosL[m];

    while(fabs(L-L0) > 1.e-8 )
     {
      if(L < L0)  //poin lies in left side
       {
        r = m;
       }
      else
       {
        l=m;
       }

      m= (l+r)/2;
      L0 = IntlPosL[m];
     } //  while ( 

    RootPtIndex[i] = m;
   }

  FeSpace_L1->SetIntlPtIndexOfPts(RootPtIndex);
  FeSpace_L1->SetN_RootNodalPts(N_RootPts);

  N_L1nodal = N_RootPts;

 //  cout << N_AllLocalPoints << " N_RootPts  "  << N_RootPts << endl;
 //   for(i=0; i<N_RootPts; i++)
 //   cout << i << " L: "  << IntlPosL[i] << endl;
 // exit(0);
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


/** transfer sol from internal to space */
void GetSolFromNodalPtVales(int N_Levels, int MaxN_PtsForNodal, TFEFunction2D **ScalarFunctions, double *Sol, int N_U, double *ValuesT, int N_Xpos)
{
  int i, ii, j, k, l, m, N_Cells, N_Points, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers, *DOF, disp;
  int PolynomialDegree, ApproxOrder, *IntlPtIndexOfPts, *PtIndexLoc;

  double *xi, *eta, *NodalValues, s, *Values_Level, *Sol_Level;
  double *val, *PtValues, maxval=0.;
  double FunctionalValues[MaxN_PointsForNodal2D];

  TFESpace2D *PBE_Spaces;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId_PBE;
  TFE2D *Element_PBE;
  TNodalFunctional2D *nf_PBE;

  // assume all  ScalarFunctions use same fespace2D
  PBE_Spaces = ScalarFunctions[0]->GetFESpace2D();
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();
  IntlPtIndexOfPts = PBE_Spaces->GetIntlPtIndexOfPts();
  val = new double[MaxN_PtsForNodal*N_Levels];

  // assume that both fespace2D and PBE_Spaces use same coll
  Coll = PBE_Spaces->GetCollection();
  N_Cells = Coll->GetN_Cells();

  disp = 0;

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    nf_PBE = Element_PBE->GetNodalFunctional2D();
    nf_PBE->GetPointsForAll(N_Points, xi, eta);
    DOF = GlobalNumbers + BeginIndex[i];
    N_LocalDOFs = Element_PBE->GetN_DOF();
    PtIndexLoc = IntlPtIndexOfPts + disp;

    // collect point values for all level
    for(j=0;j<N_Points;j++)
    {
      k = PtIndexLoc[j];
      Values_Level  =  ValuesT + k*N_Levels;

      for(ii=0;ii<N_Levels;ii++)
      {
        val[ii*N_Points + j] = Values_Level[ii];
        if(maxval<Values_Level[ii]) maxval=Values_Level[ii];
      }
    }     // for(j=0;j<N_Points;j++)

    for(ii=0;ii<N_Levels;ii++)
    {
      PtValues = val + ii*N_Points;
      Sol_Level = Sol+ ii*N_U;

      nf_PBE->GetAllFunctionals(Coll, (TGridCell *)cell, PtValues, FunctionalValues);

      for(j=0;j<N_LocalDOFs;j++)
        Sol_Level[DOF[j]] = FunctionalValues[j];
    }

    disp +=N_Points;
  }                                               // for(i=0;i<N_Cells;i++)

  delete [] val;
  //    cout << " GetSolFromQuadPtVales " << maxval << endl;
}    //void GetSolFromNodalPtVales(


void GetSolFromQuadPtVales(int N_Levels, int MaxN_PtsForNodal, TFEFunction2D **ScalarFunctions,
double *Sol, int N_U, double *ValuesT, int N_Xpos)
{
  int i, j, k, l, m, n, N_Cells, N_Points, N_LocalDOFs, N_Sets=1;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int PolynomialDegree, *N_IncidentArray;

  double *xi, *eta, *weights, *RHS, maxval=0.;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];
  double **origvalues, *sol, *org, w, val, G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];;

  bool Needs2ndDer[1];

  TFESpace2D *PBE_Spaces;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId_PBE;
  TFE2D *Element_PBE;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  TRefTrans2D *F_K;
  BaseFunct2D BaseFunct;

  // assume all  ScalarFunctions use same fespace2D
  PBE_Spaces = ScalarFunctions[0]->GetFESpace2D();
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();

  Needs2ndDer[0] = FALSE;
  RHS = new double[N_Levels*MaxN_QuadPoints_2D];
  N_IncidentArray = new int[N_U];
  memset(N_IncidentArray, 0, N_U*SizeOfInt);
  memset(Sol, 0, N_U*N_Levels*SizeOfDouble);

  // assume that both fespace2D and PBE_Spaces use same coll
  Coll = PBE_Spaces->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  m=0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);

    N_LocalDOFs = Element_PBE->GetN_DOF();
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId_PBE);
    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_PBE);
    BaseFunct = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(FEId_PBE);

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TQuadBilinear *)F_K)->GetOrigValues(BaseFunct, N_Points, xi, eta, N_LocalDOFs, QuadFormula);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TTriaAffin *)F_K)->GetOrigValues(BaseFunct, N_Points, xi, eta, N_LocalDOFs, QuadFormula);
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;
    }                                             // endswitch

    //cout << "QuadFormula: " << QuadFormula << endl;

    origvalues=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    DOF = GlobalNumbers + BeginIndex[i];

    memset(RHS, 0, N_Levels*N_LocalDOFs*SizeOfDouble);
    memset(G, 0, N_LocalDOFs*N_LocalDOFs*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      w = AbsDetjk[j]*weights[j];
      org = origvalues[j];
      sol = ValuesT+(m*N_Levels);
      m++;

      for(k=0;k<N_LocalDOFs;k++)
      {
        val = w*org[k];

        // multiple r.h.s
        for(l=0;l<N_Levels;l++)
        {
          RHS[k*N_Levels + l] += val*sol[l];
          if(maxval<sol[l]) maxval=sol[l];
        }

        //single matrix
        for(n=0;n<N_LocalDOFs;n++)
        {
          G[k*N_LocalDOFs + n] += val*org[n];
        }                                         // for(n=0;n<N_LocalDOFs;n++)
      }                                           //for(k=0;k<N_Levels;k++)
    }                                             // for(j=0;j<N_Points;j++)

    //    for(k=0;k<N_LocalDOFs;k++)
    //     for(n=0;n<N_LocalDOFs;n++)
    //     cout<< "(" << k<< " , " << n << ") = " << G[k*N_LocalDOFs + n] << " AT " <<  G[n*N_LocalDOFs + k] <<endl;

    SolveMultipleSystemsNew(G, RHS, N_LocalDOFs, N_LocalDOFs, N_Levels, N_Levels);

    for(j=0;j<N_Levels;j++)
    {
      sol = Sol+(N_U*j) ;

      for(k=0;k<N_LocalDOFs;k++)
      {
        l = DOF[k];
        if(j==0)
        {
          N_IncidentArray[l]++;
          //           OwnDof[l] = TRUE;
        }
        sol[l] += RHS[k*N_Levels + j];
      }
    }                                             // for(j=0;j<N_Levels;j++)
  }                                               // for(i=0;i<N_Cells;i++)

  for(j=0;j<N_Levels;j++)
  {
    sol = Sol+(N_U*j);

    for(k=0;k<N_U;k++)
    {
      //         if(j==0 && OwnDof[k])
      if(j==0)
      {
        if(N_IncidentArray[k]==0)
        {
          cout<< "error in GetSolFromQuadPtVales " <<endl;
          exit(0);
        }
      }
      sol[k] /= double(N_IncidentArray[k]);

    }                                             // for(k=0;k<N_LocalDOFs;k++)
  }                                               // for(j=0;j<N_Levels;j++)

  delete [] RHS;
  delete [] N_IncidentArray;

  //  cout << " GetSolFromQuadPtVales " << maxval << endl;
}   //void GetSolFromQuadPtVales(


void GetSolFromQuadPtVales(TFEFunction2D *ScalarFunction, double *Values, TFESpace2D *PBE_Spaces)
{
  int i, j, k, l, m, n, N_Cells, N_Points, N_LocalDOFs, N_U;
  int *BeginIndex, *GlobalNumbers, *DOF, PolynomialDegree, *N_IncidentArray;

  double *xi, *eta, *weights, RHS[MaxN_BaseFunctions2D], G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double val, sol, PtValues[MaxN_PointsForNodal2D], *Sol;
  double FunctionalValues[MaxN_PointsForNodal2D];
  double **origvalues, *org;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D], w;

  //   bool  *OwnDof;

  TFESpace2D *FE_Space;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId, FEId_PBE;
  TFE2D *Element, *Element_PBE;
  TNodalFunctional2D *nf;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  BaseFunct2D BaseFunct;
  TRefTrans2D *F_K;

  // assume all  ScalarFunctions use same fespace2D
  Sol =  ScalarFunction->GetValues();
  FE_Space = ScalarFunction->GetFESpace2D();
  BeginIndex = FE_Space->GetBeginIndex();
  GlobalNumbers = FE_Space->GetGlobalNumbers();
  N_U = FE_Space->GetN_DegreesOfFreedom();
  N_IncidentArray = new int[N_U];
  memset(N_IncidentArray, 0, N_U*SizeOfInt);
  memset(Sol, 0, N_U*SizeOfDouble);

  // assume that both fespace2D and FE_Space use same coll
  Coll = FE_Space->GetCollection();
  N_Cells = Coll->GetN_Cells();

  m = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId_PBE);
    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_PBE);

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;
    }                                             // endswitch

    FEId = FE_Space->GetFE2D(i, cell);
    BaseFunct = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(FEId);

    origvalues=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    DOF = GlobalNumbers + BeginIndex[i];

    memset(RHS, 0, N_LocalDOFs*SizeOfDouble);
    memset(G, 0, N_LocalDOFs*N_LocalDOFs*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      w = AbsDetjk[j]*weights[j];
      org = origvalues[j];
      sol = Values[m];
      m++;

      for(k=0;k<N_LocalDOFs;k++)
      {
        val = w*org[k];
        RHS[k] +=sol*val;

        //single matrix
        for(n=0;n<N_LocalDOFs;n++)
        {
          G[k*N_LocalDOFs + n] += val*org[n];
        }                                         // for(n=0;n<N_LocalDOFs;n++)
      }                                           // for(k=0;k<N_LocalDOFs;k++)
    }                                             // for(j=0;j<N_Points;j++)

    SolveMultipleSystemsNew(G, RHS, N_LocalDOFs, N_LocalDOFs, 1, 1);

    for(k=0;k<N_LocalDOFs;k++)
    {
      l = DOF[k];
      N_IncidentArray[l]++;
      //       OwnDof[l] = TRUE;
      Sol[l] += RHS[k];
    }
  }                                               // for(i=0;i<N_Cells;i++)

  for(k=0;k<N_U;k++)
  {
    //     if(N_IncidentArray[k]==0 && OwnDof[k])
    if(N_IncidentArray[k]==0)
    {
      cout<< "error in GetSolFromQuadPtVales " <<endl;
      exit(0);
    }

    Sol[k] /= double(N_IncidentArray[k]);
  }                                               // for(k=0;k<N_LocalDOFs;k++)


} //void GetSolFromQuadPtVales(



void GetSolFromNodalPtVales(TFEFunction2D *ScalarFunction, double *Values, TFESpace2D *PBE_Spaces)
{
  int i, ii, j, k, l, m, N_Cells, N_Points, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers, *DOF, disp, *IntlPtIndexOfPts, *PtIndexLoc;

  double *xi, *eta;
  double val[MaxN_PointsForNodal2D], PtValues[MaxN_PointsForNodal2D], *Sol;
  double FunctionalValues[MaxN_PointsForNodal2D];

  TFESpace2D *FE_Space;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TNodalFunctional2D *nf;

  // assume all  ScalarFunctions use same fespace2D
  Sol =  ScalarFunction->GetValues();
  FE_Space = ScalarFunction->GetFESpace2D();
  BeginIndex = FE_Space->GetBeginIndex();
  GlobalNumbers = FE_Space->GetGlobalNumbers();

  IntlPtIndexOfPts = PBE_Spaces->GetIntlPtIndexOfPts();

  // assume that both fespace2D and FE_Space use same coll
  Coll = FE_Space->GetCollection();
  N_Cells = Coll->GetN_Cells();

  disp = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FE_Space->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);
    DOF = GlobalNumbers + BeginIndex[i];
    N_LocalDOFs = Element->GetN_DOF();

    PtIndexLoc = IntlPtIndexOfPts+disp;

    for(j=0;j<N_LocalDOFs;j++)
    {
      k = PtIndexLoc[j];
      PtValues[j] = Values[k];
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PtValues, FunctionalValues);

    for(j=0;j<N_LocalDOFs;j++)
      Sol[DOF[j]] = FunctionalValues[j];

    disp +=N_Points;
  }                                               // for(i=0;i<N_Cells;i++)

}                                                 //void GetSolFromNodalPtVales(


void  L_Nodal2Sol(TFESpace1D *FESpace1D, double *Sol_QuadIntl, int N_XLocPoints, 
                  int N_Levels, double *Sol_AllL)
{

  int i,j,k,l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs, *IncidentArray;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, ii;
  int *DOF, *Index, *IndexOfNodalPts, disp;

  double *xi, *eta, *Values_Level, *sol;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double *PointValues, *PtVal;
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TNodalFunctional1D *nf;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  IndexOfNodalPts = FESpace1D->GetIntlPtIndexOfPts();

  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  PointValues = new double [MaxN_PointsForNodal1D*N_XLocPoints];
  IncidentArray = new int [N_DOFs];
  memset(Sol_AllL, 0, SizeOfDouble*N_DOFs*N_XLocPoints);
  memset(IncidentArray, 0, SizeOfInt*N_DOFs);

  disp = 0; 

  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    DOF = GlobalNumbers+BeginIndex[i];

    Index = IndexOfNodalPts+disp;

    for(ii=0;ii<N_XLocPoints;ii++)
     {
      Values_Level  =  Sol_QuadIntl + ii*N_Levels;

      for(j=0;j<N_Points;j++)
      {
       k = Index[j]; // corresponding L level
       PointValues[ii*N_Points + j] = Values_Level[k];
      }
     }

    for(ii=0;ii<N_XLocPoints;ii++)
     {
      PtVal = PointValues + ii*N_Points;
      nf->GetAllFunctionals(PtVal, FunctionalValues);

      sol = Sol_AllL + ii*N_DOFs;

      for(j=0;j<N_LocalDOFs;j++)
       {
        sol[DOF[j]] += FunctionalValues[j];

        if(ii==0)
         IncidentArray[DOF[j]] ++;
       }
     }

    disp +=N_Points;
  } // for(i=0;i<N_Cells

  for(ii=0;ii<N_XLocPoints;ii++)
    {
     for(i=0;i<N_DOFs;i++)
      {
       if(ii==0)
        if(IncidentArray[i] == 0)
         {
          cout << "Error in L_Nodal2Sol : "<< IncidentArray[i] << endl;
          exit(0);
         }
       Sol_AllL[ii*N_DOFs + i] /= (double)IncidentArray[i];
      } // for(i=0;i<N_DOFs;i
    } // for(ii=0;ii<N_XLocPoints;i

  delete [] PointValues;
  delete [] IncidentArray;
}


void  L_Sol2Nodal(TFESpace1D *FESpace1D, double *Sol_AllL, int N_Levels,  double *Sol_NodalPts, int N_XPoints)
{
  int i,j,k,l, ii;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, disp;
  int *DOF, *IndexArray, *NodalPtIndex;

  double *xi, *eta, *sol, *sol_Nodal, val;
  double Z[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double BasisValues[MaxN_PointsForNodal1D][MaxN_BaseFunctions1D];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();
  NodalPtIndex = FESpace1D->GetIntlPtIndexOfPts();

  IndexArray = new int[N_Levels];
  memset(IndexArray, 0, SizeOfInt*N_Levels);
  memset(Sol_NodalPts , 0, SizeOfDouble*N_XPoints*N_Levels);

  disp = 0;
  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Z, AbsDetjk);
    DOF = GlobalNumbers + BeginIndex[i];

     for(j=0;j<N_Points;j++)
       bf->GetDerivatives(D0, xi[j], BasisValues[j]);

     for(ii=0;ii<N_XPoints;ii++)
      {
       sol = Sol_AllL + ii*N_DOFs;
       sol_Nodal =  Sol_NodalPts + ii*N_Levels; 

//        for(l=0;l<N_LocalDOFs;l++)
//         cout << " val " << sol[DOF[l]] << endl;

       for(j=0;j<N_Points;j++)
        {
         k = NodalPtIndex[disp + j]; // find the right L level
         val = 0.;

          for(l=0;l<N_LocalDOFs;l++)
           val += sol[DOF[l]]*BasisValues[j][l];

         sol_Nodal[k] += val;
         if(ii==0)
          IndexArray[k]++;
        } // for(j=0;j<N_Points
      } // for(ii=0

     disp +=N_Points;
   } // for(i=0; i<N_Cells; i


   for(ii=0;ii<N_XPoints;ii++)
    {
     sol_Nodal =  Sol_NodalPts + ii*N_Levels; 
     for(i=0;i<N_Levels;i++)
      {
       if(ii==0)
        if(IndexArray[i] == 0)
         {
          cout << "Error in L_Sol2Nodal : "<< IndexArray[i] << endl;
          exit(0);
         }

       sol_Nodal[i] /= (double)IndexArray[i];
      }
    }
}



void GetOSError(int N_L1nodal, TFEFunction2D **ScalarFunctions,  TFESpace1D *FESpace1D,
double *Sol_L2nodalL1nodalXdof, int N_U, double *errors, DoubleFunctND *ExactFunct)
{
  int i,j,k,l, m, n, N_Cells, N_V, N_LDof;
  int N_DOFs, N_LocalDOFs, N_BaseFunct_Intl;
  int *BeginIndex, *GlobalNumbers, *DOF, *BeginIndex_Intl, *GlobalNumbers_Intl, *DOF_Intl;
  int N_LocalUsedElements, N_Points, *N_BaseFunct, N_Cells_Intl;
  int L, N_LinePoints, N_Sets=1;

  double *weights, *xi, *eta, *Sol_AllL;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Sol_QuadIntl, *sol;
  double **OrigFEValues, *Orig, value, *InternalError, *InternalErrorgrad;
  double *LineWeights, *zeta, Z[MaxN_QuadPoints_1D], LineAbsDetjk[MaxN_QuadPoints_1D];;
  double **origvaluesD0, **origvaluesD1, *orgD0, *orgD1, Mult;
  double x, y, z, Exact[5], *sol_QI, valuegrad, Xi[3];

  TBaseCell *cell, *Cell_Intl;
  TCollection *Coll, *Coll_Intl;
  TFESpace2D *FESpace2D;
  TFE2D *Element;
  FE2D LocalUsedElements[1], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  FE1D FEId_Intl;
  TFE1D *Element_Intl;
  TBaseFunct1D *bf_Intl;
  BaseFunct1D BaseFunct_ID_Intl, BaseFunct_Intl[1];
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans1D *F_K;

  bool *SecondDer;
  bool Needs2ndDer[1];

  FESpace2D = ScalarFunctions[0]->GetFESpace2D();
  Coll = FESpace2D->GetCollection();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_Cells = Coll->GetN_Cells();

  Coll_Intl = FESpace1D->GetCollection();
  BeginIndex_Intl = FESpace1D->GetBeginIndex();
  GlobalNumbers_Intl = FESpace1D->GetGlobalNumbers();
  N_Cells_Intl = Coll_Intl->GetN_Cells();
  N_LDof = FESpace1D->GetN_DegreesOfFreedom();
  N_LocalUsedElements = 1;
  SecondDer = new bool[1];
  SecondDer[0] = FALSE;
  Needs2ndDer[0] = FALSE;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  //first check how many quad pts  in the cell
  cell = Coll->GetCell(0);
  LocalUsedElements[0] = FESpace2D->GetFE2D(0, cell);
  Element = TFEDatabase2D::GetFE2D(LocalUsedElements[0]);
  TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, Coll, cell, SecondDer,
    N_Points, xi, eta, weights, X, Y, AbsDetjk);

  Sol_QuadIntl = new double[N_Points*N_L1nodal];
  Sol_AllL = new double[N_Points*N_LDof];
  InternalError = new double[N_Points];
  InternalErrorgrad = new double[N_Points];

  errors[0] = 0.;
  errors[1] = 0.;

  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    LocalUsedElements[0] = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(LocalUsedElements[0]);

  // ####################################################################
  // calculate values on original element
  // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);
    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace2D->GetFE2D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_LocalDOFs = N_BaseFunct[CurrentElement];
    DOF = GlobalNumbers+BeginIndex[i];

    OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);

    // find values at all quad points
    memset(Sol_QuadIntl, 0, N_L1nodal*N_Points*SizeOfDouble);
    for(j=0; j<N_L1nodal; j++)
    {
      sol = Sol_L2nodalL1nodalXdof+j*N_U;

      for(k=0; k<N_Points; k++)
      {
        Orig = OrigFEValues[k];
        value = 0.;

        for(l=0; l<N_LocalDOFs; l++)
          value += sol[DOF[l]]*Orig[l];

        Sol_QuadIntl[k*N_L1nodal  +  j] = value;
      }
    }      //for(j=0; j<N_L1nodal; j++)

   L_Nodal2Sol(FESpace1D, Sol_QuadIntl, N_Points, N_L1nodal, Sol_AllL);

    memset(InternalError, 0, N_Points*SizeOfDouble);
    memset(InternalErrorgrad, 0, N_Points*SizeOfDouble);

    for(j=0; j<N_Points; j++)
    {
      x = X[j];
      y = Y[j];
//       sol_QI = Sol_QuadIntl+j*N_L1nodal;
      sol_QI = Sol_AllL + j*N_LDof;

      for(k=0; k<N_Cells_Intl; k++)
      {
        Cell_Intl = Coll_Intl->GetCell(k);
        FEId_Intl = FESpace1D->GetFE1D(k, Cell_Intl);
        Element_Intl = TFEDatabase2D::GetFE1D(FEId_Intl);
        bf_Intl = Element_Intl->GetBaseFunct1D();
        N_BaseFunct_Intl = Element_Intl->GetN_DOF();
        BaseFunct_ID_Intl = Element_Intl->GetBaseFunct1D_ID();
        DOF_Intl = GlobalNumbers_Intl+BeginIndex_Intl[k];

        L = bf_Intl->GetPolynomialDegree();
        LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*L);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
        ((TLineAffin *)F_K)->SetCell(Cell_Intl);
        ((TLineAffin *)F_K)->GetOrigFromRef(N_LinePoints, zeta, Z, LineAbsDetjk);

        BaseFunct_Intl[0] = BaseFunct_ID_Intl;
        ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct_Intl, N_LinePoints, zeta,  LineQuadFormula,  Needs2ndDer);

        origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_Intl, D0);
        origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_Intl, D1);
        DOF_Intl = GlobalNumbers_Intl + BeginIndex_Intl[k];

        for(l=0; l<N_LinePoints; l++)
        {
          Mult = 0.5*LineWeights[l]*LineAbsDetjk[l];
          orgD0 = origvaluesD0[l];
          orgD1 = origvaluesD1[l];
          z = Z[l];

          //find sol at this point
          value = 0.; valuegrad= 0.;
          for(m=0; m<N_BaseFunct_Intl; m++)
          {
            value += sol_QI[DOF_Intl[m]]*orgD0[m];
            valuegrad  += sol_QI[DOF_Intl[m]]*orgD1[m];
          }
          Xi[0] = x;
          Xi[1] = y;
	        Xi[2] = z;
          
          ExactFunct(3,  Xi, Exact);
          InternalError[j] += Mult*(Exact[0]-value)*(Exact[0]-value);
          InternalErrorgrad[j] += Mult*(Exact[1]-valuegrad)*(Exact[1]-valuegrad);

//           if(j==0 && l==0)
//            cout <<  j << " " << k << " " << Exact[0] << " " << value << endl;

        }  //for(l=0; l
      } //for(k=0;

      Mult = weights[j]*AbsDetjk[j];
      errors[0] += Mult*InternalError[j];
      errors[1] += Mult*InternalErrorgrad[j];;  
    }     // for(j=0;
  }      //for(i=0; i<N_Cells; i++)

  delete [] Sol_QuadIntl;
  delete [] InternalError;

}  //void GetOSError





void AssembleFD_M_StartPtDirichlet(TSquareMatrix1D *M)
{
  int i, j, N, N_Values; 
  int begin, end, *Row, *KCol;
  double *ValuesM;
  
  
  N = M->GetN_Columns();
  Row = M->GetRowPtr();
  KCol = M->GetKCol();
  ValuesM = M->GetEntries();  
   
  M->Reset();
    
    for(i=0;i<N;i++)
    {
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
       if(i==KCol[j])
        ValuesM[j] = 1.;
      }   // endfor j
    }    
  
}


void AssembleM_StartPtDirichlet(TFESpace1D *FeSpace, TSquareMatrix1D *M)
{
  int i, j, k, l, N_Cells, N_BaseFunct, N_U, dGDisc;
  int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
  int TestDOF, begin, end, *RowPtr, *KCol;

  double *Weights, *zeta, X[20], AbsDetjk[20];
  double LocMatrixM[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
  double **origvaluesD0, **origvaluesD1, Mult;
  double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesM;
  double x=0., len =0.;

  bool Needs2ndDer[1];

  TBaseCell *Cell;
  FE1D FEId;
  TFE1D *Element;
  TBaseFunct1D *bf;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans1D *F_K;
  BaseFunct1D BaseFunct_ID, BaseFunct[1];
  TCollection *Coll;
  BoundCond BDType;
  BoundCond cond_Lmin, cond_Lmax;

  BoundCondition_LminLMax_L1(cond_Lmin, cond_Lmax);

  dGDisc=FeSpace->IsDGSpace();
  Coll = FeSpace->GetCollection();
  N_U = FeSpace->GetN_DegreesOfFreedom();
  GlobalNumbers = FeSpace->GetGlobalNumbers();
  BeginIndex = FeSpace->GetBeginIndex();
  RowPtr = M->GetRowPtr();
  KCol = M->GetKCol();
  ValuesM = M->GetEntries();

  // all QuadPts in a cell use same FEspace in internal direction
  N_Cells = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  for(i=0; i<N_Cells; i++)
  {
    Cell = Coll->GetCell(i);
    FEId = FeSpace->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();

    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, X, AbsDetjk);

    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);

    memset(LocMatrixM, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
    {
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];

      //cout<< " zeta[j] " << zeta[j]  <<endl;
      //len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
      {
        test0  = orgD0[k];
        //cout<< " uref " << test0  <<endl;
        for(l=0;l<N_BaseFunct;l++)
        {
          ansatz0  = orgD0[l];
          LocMatrixM[k*N_BaseFunct + l] += (Mult*ansatz0*test0);
        }
      }
    }

    //   add to global matrices
    for(j=0;j<N_BaseFunct;j++)
    {
      TestDOF = DOF[j];

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
      {
        for(l=0;l<N_BaseFunct;l++)
        {
          if(KCol[k] == DOF[l])
          {
            ValuesM[k] +=LocMatrixM[j*N_BaseFunct + l];
            break;
          }
        }                                         // for(m=0;m<N_BaseFunct_low
      }                                           // for(n=begin;n<end;n++)
    }                                             // for(l=0;l<N_BaseFunct_low
  }                                               // for(i=0; i<N_Cells; i++)

  //update boundary data
  // starting point:  in PBS
  // end point: in PBS

  if(cond_Lmin==DIRICHLET && !dGDisc)
  {
    begin = RowPtr[0];
    end = RowPtr[1];

    for(k=begin;k<end;k++)
    {
      if(KCol[k] == 0 )
        { ValuesM[k] = 1.; }
        else
          { ValuesM[k] = 0.; }
    }
  }                                             //if(cond_Lmin==DIRICHLET)

  if(cond_Lmax==DIRICHLET && !dGDisc)
  {
    begin = RowPtr[N_U-1];
    end = RowPtr[N_U];

    for(k=begin;k<end;k++)
    {
      if(KCol[k] == N_U-1 )
        { ValuesM[k] = 1.; }
        else
          { ValuesM[k] = 0.; }
    }
  }                                             //if(cond_Lmin==DIRICHLET)

  // cout<< " len " << len << endl;
  // exit(0);
  // //print matrix
  //   for(j=0;j<N_U;j++)
  //    {
  //     begin = RowPtr[j];
  //     end = RowPtr[j+1];
  //     for(k=begin;k<end;k++)
  //      {
  //       cout << "M(" << j << ", "<< KCol[k] << ") = " << ValuesM[k] <<endl;
  //      }
  //     cout<<endl;
  //    }
  //    exit(0);
}


 /** transfer sol from spatial to internal */
void GetPtsValuesForIntl(int N_Levels, TFESpace2D *PBE_Spaces, double *Sol, int N_U, double *ValuesT, int N_Xpos)
{
  int i, ii, j, k, l, m, N_Cells, N_Points, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int PolynomialDegree, ApproxOrder, *IntlPtIndexOfPts, *PtIndexLoc, disp;
  int *IncidentArray, *Incident;

  double *xi, *eta, *NodalValues, s, *Values_Level;
  double BasisValues[MaxN_BaseFunctions2D], maxval=0.;

  // TFESpace2D *PBE_Spaces;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId_PBE;
  TFE2D *Element_PBE;
  TBaseFunct2D *bf_PBE;
  TNodalFunctional2D *nf_PBE;

  // assume all  ScalarFunctions use same fespace2D
  // PBE_Spaces = ScalarFunctions[0]->GetFESpace2D();
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();

  IntlPtIndexOfPts = PBE_Spaces->GetIntlPtIndexOfPts();
  IncidentArray = new int[N_Levels*N_Xpos];

  memset(ValuesT, 0, N_Levels*N_Xpos*SizeOfDouble);
  memset(IncidentArray, 0, N_Levels*N_Xpos*SizeOfInt);

  // assume that both fespace2D and PBE_Spaces use same coll
  Coll = PBE_Spaces->GetCollection();
  N_Cells = Coll->GetN_Cells();

  disp = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    nf_PBE = Element_PBE->GetNodalFunctional2D();
    nf_PBE->GetPointsForAll(N_Points, xi, eta);
    bf_PBE = Element_PBE->GetBaseFunct2D();
    N_LocalDOFs = Element_PBE->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];
    PtIndexLoc = IntlPtIndexOfPts + disp;

    for(j=0;j<N_Points;j++)
    {
      bf_PBE->GetDerivatives(D00, xi[j], eta[j], BasisValues);

      k = PtIndexLoc[j];

      Incident = IncidentArray + k*N_Levels;
      Values_Level  =  ValuesT + k*N_Levels;

      for(ii=0;ii<N_Levels;ii++)
      {
        Incident[ii]++;
        NodalValues =  Sol + ii*N_U;
        for(l=0;l<N_LocalDOFs;l++)
        {
          m = DOF[l];
          s = NodalValues[m];
          Values_Level[ii] += BasisValues[l]*s;
        }
      }                                           // for(ii=0;ii<N_Levels;ii++)
    }                                             // for(j=0;j<N_Points;j++)

    disp +=N_Points;
  }                                               // for(i=0;i<N_Cells;i++)

  for(i=0;i<N_Xpos;i++)
  {
    Incident = IncidentArray + i*N_Levels;
    Values_Level  =  ValuesT + i*N_Levels;

    for(j=0;j<N_Levels;j++)
    {
      Values_Level[j] /= (double)Incident[j];
      // if(maxval < Values_Level[j]) maxval = Values_Level[j];
    }
  }

  delete [] IncidentArray;

  // cout <<  " GetPtsValuesForIntl_Nodal " << maxval <<endl;
  //exit(0);
}

void GetPtsValuesForIntl_Quad(int N_Levels, TFESpace2D *PBE_Spaces, double *Sol, int N_U, double *ValuesT, int N_Xpos)
{
  int i, j, k, l, m, N_Cells, N_Points, N_LocalDOFs, N_Sets=1;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int PolynomialDegree;

  double *xi, *eta, *weights, maxval=0.;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];
  double **origvalues, *sol, *org, val;

  bool Needs2ndDer[1];

  // TFESpace2D *PBE_Spaces;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId_PBE;
  TFE2D *Element_PBE;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  TRefTrans2D *F_K;
  BaseFunct2D BaseFunct;

  // assume all  ScalarFunctions use same fespace2D
  // PBE_Spaces = ScalarFunctions[0]->GetFESpace2D();
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();

  Needs2ndDer[0] = FALSE;

  // assume that both fespace2D and PBE_Spaces use same coll
  Coll = PBE_Spaces->GetCollection();
  N_Cells = Coll->GetN_Cells();
  m = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    N_LocalDOFs = Element_PBE->GetN_DOF();
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId_PBE);
    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_PBE);
    BaseFunct = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(FEId_PBE);

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TQuadBilinear *)F_K)->GetOrigValues(BaseFunct, N_Points, xi, eta, N_LocalDOFs, QuadFormula);
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TTriaAffin *)F_K)->GetOrigValues(BaseFunct, N_Points, xi, eta, N_LocalDOFs, QuadFormula);
        break;
    }                                             // endswitch

    //cout << "QuadFormula: " << QuadFormula << endl;

    origvalues=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0; j<N_Levels; j++)
    {
      // sol = ScalarFunctions[j]->GetValues();
      sol =  Sol + j*N_U;
      for(k=0; k<N_Points; k++)
      {
        org = origvalues[k];
        val = 0.;
        for(l=0;l<N_LocalDOFs;l++)
        {
          val += org[l]*sol[ DOF[l] ];            // for M
        }

        ValuesT[(m+k)*N_Levels + j] = val;

        if(maxval<val) maxval=val;

      }                                           //  for(k=0; k<N_Points; k++)
    }                                             // for(j=0; j<N_Levels; j++)

    m += N_Points;
  }                                               //  for(i=0;i<N_Cells;i++)

  // cout <<  " GetPtsValuesForIntl_Quad " << maxval <<endl;

  //    cout <<  m << " N_Xpos " << N_Xpos <<endl;
  // exit(0);
}


void GetSol_L2_X_L1_and_Sol_X_IntOmegaL1L2(TFESpace1D *FeSpace_L1, TFESpace1D *FeSpace_L2, TFESpace2D *FeSpace_X, 
                            int N_Xpos, int N_L1nodal, int N_L2nodal, double *Sol_L2nodalL1nodalXdof, 
                            double *Sol_L2nodalXposL1nodal, double *Sol_L2nodalXposL1dof, double *Sol_Xpos_IntOmegaL2L1,
                            double *IntValue)
{
 double *Sol_XposL1nodal, *Sol_L1nodalXdof, *Sol_XposL1dof, *Sol_XposL1quadL2Dof;
 int N_DOF = FeSpace_X->GetN_DegreesOfFreedom();
 int N_L1DOF = FeSpace_L1->GetN_DegreesOfFreedom();
 int N_L2DOF = FeSpace_L2->GetN_DegreesOfFreedom();
 int i, ii, j, k;
 int i2, ii2, j2, k2;
 FE1D LocalUsedElements[1], CurrentElement;

  TCollection *L1Coll = FeSpace_L1->GetCollection();
  int *L1BeginIndex = FeSpace_L1->GetBeginIndex();
  int *L1GlobalNumbers = FeSpace_L1->GetGlobalNumbers();
  int N_L1Cells = L1Coll->GetN_Cells();
  TCollection *L2Coll = FeSpace_L2->GetCollection();
  int *L2BeginIndex = FeSpace_L2->GetBeginIndex();
  int *L2GlobalNumbers = FeSpace_L2->GetGlobalNumbers();
  int N_L2Cells = L2Coll->GetN_Cells();

  bool Needs2ndDer[1], SecondDer[1];
  SecondDer[0] = FALSE;
  Needs2ndDer[0] = FALSE;

 for(int ii=0; ii<N_L2nodal; ii++)
  {
   Sol_XposL1dof = Sol_L2nodalXposL1dof+(ii*N_L1DOF*N_Xpos);

   #if defined(__SPATIAL__) || defined(__L2DIR__)
   /** transfer sol from spatial to internal */   
   Sol_L1nodalXdof = Sol_L2nodalL1nodalXdof+(ii*N_L1nodal*N_DOF);       
   Sol_XposL1nodal = Sol_L2nodalXposL1nodal+(ii*N_L1nodal*N_Xpos);

   /** transfer sol from spatial to internal */
   #ifdef __OSNODALPT__
   GetPtsValuesForIntl(N_L1nodal, FeSpace_X, Sol_L1nodalXdof, N_DOF, Sol_XposL1nodal, N_Xpos);
   #else
   GetPtsValuesForIntl_Quad(N_L1nodal, FeSpace_X,  Sol_L1nodalXdof, N_DOF, Sol_XposL1nodal, N_Xpos);
   #endif

   /** evaluate dof values from nodal points */
   L_Nodal2Sol(FeSpace_L1, Sol_XposL1nodal, N_Xpos, N_L1nodal, Sol_XposL1dof);
   #endif

  //  for(j=0; j<N_Xpos; j++)
  //   {
  //    for(k=0; k<N_L1DOF; k++)
  //     cout<< ii << "," << j << "," << k << ": Sol_XposL1dof " << Sol_XposL1dof[j*N_L1DOF + k] <<endl;
  //    cout<<endl;   
  //   }
   } // for(int ii=0; 
// exit(0);

  //first check how many quad pts  in the l1cell
  TBaseCell *cell_L2, *cell = L1Coll->GetCell(0);
  LocalUsedElements[0] = FeSpace_L1->GetFE1D(0, cell);
  TFE1D *Element = TFEDatabase2D::GetFE1D(LocalUsedElements[0]);
  TBaseFunct1D *bf = Element->GetBaseFunct1D();
  int l = bf->GetPolynomialDegree();
  QuadFormula1D LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
  TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  int N_QuadPts, N_QuadPts_L2, *L1DOF, *L2DOF;
  double *Weights, *zeta, X[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D], *Weights_L2, *zeta_L2, AbsDetjk_L2[MaxN_QuadPoints_1D];
  double Mult, value, *orgD0, **origvaluesD0;
  qf1->GetFormulaData(N_QuadPts, Weights, zeta);

  int N_BaseFunct, N_Sets=1;
  BaseFunct1D BaseFunct_ID, BaseFunct[1];
  TRefTrans1D *F_K;

  /** storage for N_QuadPts in one L1 cell */ 
  double *Sol_L1dof, *Sol_XposL1quadL2nodal = new double[N_Xpos*N_QuadPts*N_L2nodal];
  Sol_XposL1quadL2Dof = new double[N_Xpos*N_QuadPts*N_L2DOF];
  double *Sol_L1quadL2Dof, *Sol_L1quadL2nodal, *Sol_L2;
  double *IntValue_Xpos, *IntValue_XposL1quad = new double[N_Xpos*N_QuadPts]; 

  memset(IntValue, 0, N_Xpos*SizeOfDouble);
  memset(Sol_XposL1quadL2nodal, -1, N_Xpos*N_QuadPts*N_L2nodal*SizeOfDouble);

  //compute sol at L1 quadrature values
  for(i=0; i<N_L1Cells; ++i)
   {
    cell = L1Coll->GetCell(i);
    LocalUsedElements[0] = FeSpace_L1->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(LocalUsedElements[0]);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_QuadPts, Weights, zeta);
    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_QuadPts, zeta, X, AbsDetjk);
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_QuadPts, zeta,  LineQuadFormula,  Needs2ndDer);
    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    L1DOF = L1GlobalNumbers + L1BeginIndex[i];
    memset(IntValue_XposL1quad, 0, N_Xpos*N_QuadPts*SizeOfDouble);

    for(ii=0; ii<N_L2nodal; ii++)
    {
     Sol_XposL1dof = Sol_L2nodalXposL1dof+(ii*N_Xpos*N_L1DOF);

     for(j=0; j<N_Xpos; j++)
      {
       Sol_L1dof = Sol_XposL1dof +j*N_L1DOF;
       for(k=0; k<N_QuadPts; k++)
       {
        orgD0 = origvaluesD0[k];
        value = 0.;
        for(l=0;l<N_BaseFunct;l++)
        { 
         value += Sol_L1dof[L1DOF[l]]*orgD0[l]; 
        } // l
        Sol_XposL1quadL2nodal[j*N_QuadPts*N_L2nodal + k*N_L2nodal + ii] = value;
       }  // for k
     }// for j
    } // for ii

    //evaluate l2-nodal pt values to dof nodal values
    for(j=0; j<N_Xpos; j++)
     {
      Sol_L1quadL2nodal = Sol_XposL1quadL2nodal+(j*N_QuadPts*N_L2nodal);
      L_Nodal2Sol(FeSpace_L2, Sol_L1quadL2nodal, N_QuadPts, N_L2nodal, Sol_XposL1quadL2Dof+(j*N_QuadPts*N_L2DOF) );

          //  for(ii=0; ii<N_QuadPts; ii++)
          //  {
          //   for(k=0; k<N_L2nodal; k++)
          //    cout<< ii << "," <<  k << " Sol_L1quadL2nodal " << Sol_L1quadL2nodal[ii*N_L2nodal + k] <<endl;
          //    cout<<endl;  
          //  }
          //  for(ii=0; ii<N_QuadPts; ii++)
          //  {
          //   for(k=0; k<N_L2DOF; k++)
          //    cout<< ii << "," <<  k << " Sol_XposL1quadL2Dof " << Sol_XposL1quadL2Dof[ii*N_L2DOF + k] <<endl;
          //    cout<<endl;  
          //  }   
    }
    // exit(0);

   for(i2=0; i2<N_L2Cells; ++i2)
    {
     cell_L2 = L2Coll->GetCell(i2);
     LocalUsedElements[0] = FeSpace_L2->GetFE1D(i2, cell_L2);
     Element = TFEDatabase2D::GetFE1D(LocalUsedElements[0]);
     bf = Element->GetBaseFunct1D();
     N_BaseFunct = Element->GetN_DOF();
     BaseFunct_ID = Element->GetBaseFunct1D_ID();
     l = bf->GetPolynomialDegree();
     LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
     qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
     qf1->GetFormulaData(N_QuadPts_L2, Weights_L2, zeta_L2);
    //  F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
     ((TLineAffin *)F_K)->SetCell(cell_L2);
     ((TLineAffin *)F_K)->GetOrigFromRef(N_QuadPts_L2, zeta_L2, X, AbsDetjk_L2);
     BaseFunct[0] = BaseFunct_ID;
     ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_QuadPts_L2, zeta_L2,  LineQuadFormula,  Needs2ndDer);
     origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
     L2DOF = L2GlobalNumbers + L2BeginIndex[i2];
     
      for(j=0; j<N_Xpos; j++)
       {
        Sol_L1quadL2Dof = Sol_XposL1quadL2Dof+(j*N_QuadPts*N_L2DOF);        
        IntValue_Xpos = IntValue_XposL1quad+j*N_QuadPts;

        for(k=0; k<N_QuadPts; k++)
         {
          Sol_L2 = Sol_L1quadL2Dof + (k*N_L2DOF);
          value = 0.;

           // l2 integral
          for(k2=0; k2<N_QuadPts_L2; k2++)
           {
            orgD0 = origvaluesD0[k2];
            for(l=0;l<N_BaseFunct;l++)
             { 
              value += Sol_L2[L2DOF[l]]*orgD0[l]; 
            //  value += orgD0[l];              
             } // l

            value *= Weights_L2[k2]*AbsDetjk_L2[k2];
          }  // for k

          //i2 cell contribution to 
          IntValue_Xpos[k] += value;
         } // for k 
     } // for j
    } // for i2

//     for(j=0; j<N_Xpos; j++)
//      for(k=0; k<N_QuadPts; k++)
//       cout<< k << " IntValue_XposL1quad " << IntValue_XposL1quad[k] <<endl;
//     cout<<endl;
// exit(0);

    // l1 integral
    for(j=0; j<N_Xpos; j++)
     {
       IntValue_Xpos = IntValue_XposL1quad+j*N_QuadPts;
       value = 0.;
       for(k=0; k<N_QuadPts; k++)
       {
        value += IntValue_Xpos[k]*Weights[k]*AbsDetjk[k];
      }
     IntValue[j] += value;
    }
  } // for(i=0; i<N_L1Cel


   delete [] Sol_XposL1quadL2nodal;
   delete [] Sol_XposL1quadL2Dof;
   delete [] IntValue_XposL1quad;

      //  for(j=0; j<N_Xpos; j++)
    //  cout<< j << " IntValue " << IntValue[j] <<endl;
 
// exit(0);

  // cout<< " GetSol_L2_X_L1_and_Sol_X_IntOmegaL1L2 done " <<endl;
} // GetSol_L2_X_L1_and_Sol_X_IntOmegaL1L2


