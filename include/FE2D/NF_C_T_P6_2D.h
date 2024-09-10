/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_C_T_P6_2D_Xi[] = { 0, 0.16666666666666666666667, 0.33333333333333333333333, 0.5, 0.66666666666666666666667, 0.83333333333333333333333, 1, 
                                    0, 0.16666666666666666666667, 0.33333333333333333333333, 0.5, 0.66666666666666666666667, 0.83333333333333333333333,
                                    0, 0.16666666666666666666667, 0.33333333333333333333333, 0.5, 0.66666666666666666666667,
                                    0, 0.16666666666666666666667, 0.33333333333333333333333, 0.5,
                                    0, 0.16666666666666666666667, 0.33333333333333333333333,
                                    0, 0.16666666666666666666667,
                                    0 };

static double NF_C_T_P6_2D_Eta[] = {
        0,                         0,                         0,                         0,                         0,                         0,                         0,
        0.16666666666666666666667, 0.16666666666666666666667, 0.16666666666666666666667, 0.16666666666666666666667, 0.16666666666666666666667, 0.16666666666666666666667,
        0.33333333333333333333333, 0.33333333333333333333333, 0.33333333333333333333333, 0.33333333333333333333333, 0.33333333333333333333333,
        0.5,                       0.5,                       0.5,                       0.5,
        0.66666666666666666666667, 0.66666666666666666666667, 0.66666666666666666666667,
        0.83333333333333333333333, 0.83333333333333333333333,
        1 };

static double NF_C_T_P6_2D_T[] = { -1, -0.66666666666666666666667,
        -0.33333333333333333333333, 0, 0.33333333333333333333333,
        0.66666666666666666666667, 1 };

void NF_C_T_P6_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
  Functionals[6] = PointValues[6];
  Functionals[7] = PointValues[7];
  Functionals[8] = PointValues[8];
  Functionals[9] = PointValues[9];
  Functionals[10] = PointValues[10];
  Functionals[11] = PointValues[11];
  Functionals[12] = PointValues[12];
  Functionals[13] = PointValues[13];
  Functionals[14] = PointValues[14];
  Functionals[15] = PointValues[15];
  Functionals[16] = PointValues[16];
  Functionals[17] = PointValues[17];
  Functionals[18] = PointValues[18];
  Functionals[19] = PointValues[19];
  Functionals[20] = PointValues[20];
  Functionals[21] = PointValues[21];
  Functionals[22] = PointValues[22];
  Functionals[23] = PointValues[23];
  Functionals[24] = PointValues[24];
  Functionals[25] = PointValues[25];
  Functionals[26] = PointValues[26];
  Functionals[27] = PointValues[27];
}

void NF_C_T_P6_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
  Functionals[6] = PointValues[6];
}

TNodalFunctional2D *NF_C_T_P6_2D_Obj = new TNodalFunctional2D
        (NF_C_T_P6_2D, 28, 7, 28, 7, NF_C_T_P6_2D_Xi, NF_C_T_P6_2D_Eta,
         NF_C_T_P6_2D_T, NF_C_T_P6_2D_EvalAll, NF_C_T_P6_2D_EvalEdge);
