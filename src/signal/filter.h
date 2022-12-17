#ifndef _SIGNAL_FILTER_H_
#define _SIGNAL_FILTER_H_

#include <complex.h>
#include <math.h>
#include <stdlib.h>

double *signal_filter(int n_filt, double *filt, int n_a, double *a, int n_x, double *x);

#endif // _SIGNAL_FILTER_H_