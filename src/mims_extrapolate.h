#ifndef _EXTRAPOLATE_H_
#define _EXTRAPOLATE_H_

#include <stdlib.h>

#include "mims_unit.h"
#include "stats/pgamma.h"
#include "stats/spline.h"
#include "stats/smspline.h"
#include "stats/approx.h"

#ifndef max
#define max(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#endif

#ifndef min
#define min(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

typedef struct edges
{
    uint32_t n_left;
    uint32_t n_right;
    int32_t *left_start;
    int32_t *left_end;
    int32_t *right_start;
    int32_t *right_end;
} edges_t;

typedef struct values_dataframe
{
    uint32_t size;
    double *timestamps;
    double *values;
} values_dataframe_t;

dataframe_t extrapolate(dataframe_t *df, int8_t r_low, int8_t r_high, float noise_level, float k, float spar);

#endif // _EXTRAPOLATE_H_