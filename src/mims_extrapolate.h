#ifndef _EXTRAPOLATE_H_
#define _EXTRAPOLATE_H_

#include <stdlib.h>

#include "mims_unit.h"
#include "stats/pgamma.h"
#include "stats/spline.h"
#include "stats/smspline.h"
#include "stats/approx.h"

#define max(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

typedef struct edges
{
    int n_left;
    int n_right;
    int *left_start;
    int *left_end;
    int *right_start;
    int *right_end;
} edges_t;

typedef struct values_dataframe
{
    int size;
    double *timestamps;
    double *values;
} values_dataframe_t;

dataframe_t extrapolate(dataframe_t *df, double r_low, double r_high, double noise_level, double k, double spar);

#endif // _EXTRAPOLATE_H_