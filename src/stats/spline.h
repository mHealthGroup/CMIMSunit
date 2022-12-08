#ifndef _STATS_SPLINE_H_
#define _STATS_SPLINE_H_

#include <stdlib.h>

typedef struct
{
    int method;
    int n;
    double *x;
    double *y;
    double *b;
    double *c;
    double *d;
} Z_struct_t;

// Z_struct_t SplineCoef(int method, int n, double *x, int m, double *y);
// double *SplineEval(int nu, double *xout, Z_struct_t z);
double *Spline(int x_len, double *x, int y_len, double *y,
               int xout_len, double *xout, int method);

#endif // _STATS_SPLINE_H_