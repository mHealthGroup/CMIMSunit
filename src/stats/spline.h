#ifndef _STATS_SPLINE_H_
#define _STATS_SPLINE_H_

#include <stdlib.h>

double *Spline(const int x_len, double *x, const int y_len, double *y,
               const int xout_len, double *xout, const int method);

#endif // _STATS_SPLINE_H_