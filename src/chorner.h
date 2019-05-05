/*
  chorner.h

  J.J. Green 2016
*/

#ifndef CHORNER_H
#define CHORNER_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdlib.h>
#include <complex.h>
#include <math.h>

extern double crhorner(const double *p, size_t n, double x);
extern _Complex double cchorner(const _Complex double *p, size_t n,
				_Complex double x);
#ifdef __cplusplus
}
#endif

#endif
