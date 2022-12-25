#ifndef _STATS_SPLINE_H_
#define _STATS_SPLINE_H_

#include <stdlib.h>

double *Spline(int x_len, double *x, int y_len, double *y,
               int xout_len, double *xout, int method);

#endif // _STATS_SPLINE_H_