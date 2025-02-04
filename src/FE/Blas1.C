// =======================================================================
// Blas1.C
//
// Purpose:     basic routines for linear algebra
//
// Author:      Gunar Matthies          17.01.2003
//
// =======================================================================

#include <string.h>
#include <math.h>

/** return inner product (x,y) */
double Ddot(int n, double *x, double *y)
{
  register double r;
  register int i;
  register double *a, *b;

  r = 0.0;
  a = x;
  b = y;
  for(i=0;i<n;i++)
  {
    r += *a * *b;
    a++;
    b++;
  }

  return r;
}

/** return L2 Norm */
double L2Norm(int n, double *x)
{
    // Input validation
    if (n <= 0 || x == NULL)
    {
        return -1.0;
    }

    // Calculate the inner product (x,x) using Ddot
    double norm_squared = Ddot(n, x, x);

    // Handle numerical stability for very small values
    if (norm_squared < 0.0)
    {
        return 0.0;  // Return 0 for numerically unstable cases
    }

    // Return the square root of the inner product
    return sqrt(norm_squared);
}

/** return L-infinity Norm */
double LinfNorm(int n, double *x)
{
  // Loop over all elements and find the absolute maximum
  double max = -INFINITY;
  for (int i = 0; i < n; i++)
  {
    double abs_val = fabs(x[i]);
    if (abs_val > max)
    {
      max = abs_val;
    }
  }

  return max;
}

/** y := alpha*x + y */
void Daxpy(int n, double alpha, double *x, double *y)
{
  register int i;
  register double *a, *b;
  register double scal;

  a = x;
  b = y;
  scal = alpha;
  for(i=0;i<n;i++)
  {
    *b += scal * *a;
    a++;
    b++;
  }

}

/** z := alpha*x + beta*y */
void Dsum(int n, double alpha, double beta, double *x, double *y, double *z)
{
  register int i;
  register double *a, *b, *c;
  register double scal1, scal2;

  a = x;
  b = y;
  c = z;
  scal1 = alpha;
  scal2 = beta;
  for(i=0;i<n;i++)
  {
    *c = scal1 * *a + scal2 * *b;
    a++;
    b++;
    c++;
  }

}

/** z += alpha*(x + y) */
void DsumP(int n, double alpha, double *x, double *y, double *z)
{
  register int i;
  register double *a, *b, *c;
  register double scal1, scal2;

  a = x;
  b = y;
  c = z;
  scal1 = alpha;
  for(i=0;i<n;i++)
  {
    *c += scal1 * (*a + *b);
    a++;
    b++;
    c++;
  }

}

/** z += alpha*(x - y) */
void DsumS(int n, double alpha, double *x, double *y, double *z)
{
  register int i;
  register double *a, *b, *c;
  register double scal1, scal2;

  a = x;
  b = y;
  c = z;
  scal1 = alpha;
  for(i=0;i<n;i++)
  {
    *c += scal1 * (*a - *b);
    a++;
    b++;
    c++;
  }

}


/** b := a */
void Dcopy(int n, double *a, double *b)
{
  memcpy(b, a, n*sizeof(double));
}

/** x := alpha*x */
void Dscal(int n, double alpha, double *x)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a *= scal;
    a++;
  }
}

/** x := alpha+x */
void Dadd(int n, double alpha, double *x)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a += scal;
    a++;
  }
}


/** return Euclidian norm of x */
double Dnorm(int n, double *x)
{
  register double r;
  register int i;
  register double *a;
  register double z;

  a = x;
  r = 0.0;
  for(i=0; i<n; i++)
  {
    z = *a; 
    r += z*z;
    a++;
  }

  return sqrt(r);
}


