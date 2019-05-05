/*
  chorner.c

  J.J. Green 2016
*/

#include "chorner.h"

// #ifndef __STDC_IEC_559__
// #error "requires IEEE 754 floating point"
// #endif

static void two_sum(double a, double b,
		    double *ps, double *pe)
{
  double
    s = a + b,
    z = s - a,
    e = (a - (s - z)) + (b - z);

  *ps = s;
  *pe = e;
}

static void two_sum_complex(double complex a, double complex b,
			    double complex *ps, double complex *pe)
{
  double rs, re, is, ie;

  two_sum(creal(a), creal(b), &rs, &re);
  two_sum(cimag(a), cimag(b), &is, &ie);

  *ps = rs + I * is;
  *pe = re + I * ie;
}

#ifdef _WIN32
static void split(double a, double *px, double *py)
{
  double
    z = a * 0x8000001,
    x = z - (z - a),
    y = a - x;

  *px = x;
  *py = y;
}

static void two_product(double a, double b,
			double *pp, double *pe)
{
  double
    p = a * b,
    ah, al, bh, bl;

  split(a, &ah, &al);
  split(b, &bh, &bl);

  double
    e = al * bl - (((p - ah * bh) - al * bh) - ah * bl);

  *pp = p;
  *pe = e;
}
#else
static void two_product(double a, double b,
     double *pp, double *pe)
{
  *pp = a * b;
  *pe = fma(a, b, -*pp);
}
#endif

static void two_product_complex(double complex a, double complex b,
				double complex *pp,
				double complex *pe1,
				double complex *pe2,
				double complex *pe3)
{
  double
    z1, z2, z3, z4, z5, z6,
    h1, h2, h3, h4, h5, h6;

  two_product(creal(a), creal(b), &z1, &h1);
  two_product(cimag(a), cimag(b), &z2, &h2);
  two_product(creal(a), cimag(b), &z3, &h3);
  two_product(cimag(a), creal(b), &z4, &h4);
  two_sum(z1, -z2, &z5, &h5);
  two_sum(z3, z4, &z6, &h6);

  *pp  =  z5 + I * z6;
  *pe1 =  h1 + I * h3;
  *pe2 = -h2 + I * h4;
  *pe3 =  h5 + I * h6;
}

static void eft_horner(const double *p, size_t n,
		       double x,
		       double *ph,
		       double *e[static 2])
{
  double
    h = p[n-1];

  for (size_t i = 0 ; i < n - 1 ; i++)
    {
      size_t j = n - i - 2;
      double t;

      two_product(h, x, &t, e[0] + j);
      two_sum(t, p[j], &h, e[1] + j);
    }

  *ph = h;
}

static void eft_horner_complex(const double complex *p, size_t n,
			       double complex x,
			       double complex *ph,
			       double complex *e[static 4])
{
  double complex
    h = p[n-1];

  for (size_t i = 0 ; i < n - 1 ; i++)
    {
      size_t j = n - i - 2;
      double complex t;

      two_product_complex(h, x, &t, e[0] + j, e[1] + j, e[2] + j);
      two_sum_complex(t, p[j], &h, e[3] + j);
    }

  *ph = h;
}

/*
  the first argument here is "const" but mutidimensional arrays
  passes as const are treated as "pointers to const pointers to
  double" rather than "pointers to pointers to const double" so
  and we get compiler warnings which need a cast to silence ...
*/

static void horner_2sum(double *e[static 2], size_t n, double x,
			double *pr)
{
  double
    r = e[0][n-1] + e[1][n-1];

  for (size_t i = 0 ; i < n - 1 ; i++)
    {
      size_t j = n - i - 2;

      r = r * x + (e[0][j] + e[1][j]);
    }

  *pr = r;
}

// FIXME - need accurate add in the e-sum

static void horner_4sum(double complex *e[static 4], size_t n,
			double complex x, double complex *pr)
{
  double complex
    r = e[0][n-1] + e[1][n-1] + e[2][n-1] + e[3][n-1];

  for (size_t i = 0 ; i < n - 1 ; i++)
    {
      size_t j = n - i - 2;

      r = r * x + (e[0][j] + e[1][j] + e[2][j] + e[3][j]);
    }

  *pr = r;
}

static double crhorner_nontrivial(const double *p, size_t n, double x)
{
  double
    h, c,
    e0[n-1], e1[n-1], *e[2] = {e0, e1};

  eft_horner(p, n, x, &h, e);
  horner_2sum(e, n-1, x, &c);

  return h + c;
}

static double complex cchorner_nontrivial(const double complex *p, size_t n,
					  double complex x)
{
  double complex
    h, c,
    e0[n-1], e1[n-1], e2[n-1], e3[n-1],
    *e[4] = {e0, e1, e2, e3};

  eft_horner_complex(p, n, x, &h, e);
  horner_4sum(e, n-1, x, &c);

  return h + c;
}

extern double crhorner(const double *p, size_t n, double x)
{
  if (n > 1)
    return crhorner_nontrivial(p, n, x);

  if (n > 0)
    return p[0];

  return 0;
}

extern double complex cchorner(const double complex *p, size_t n,
			       double complex x)
{
  if (n > 1)
    return cchorner_nontrivial(p, n, x);

  if (n > 0)
    return p[0];

  return 0;
}
