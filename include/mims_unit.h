#ifndef _MIMS_UNIT_H_
#define _MIMS_UNIT_H_

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct
{
    int size;
    double *timestamps; // dataframe "index", seconds w/ decimal milliseconds
    int n_segments;
    int *segments;
    double *x;
    double *y;
    double *z;
    double *mims_data;
} dataframe_t;

typedef enum time_unit
{
    second,
    minute,
    hour,
    day
} time_unit_t;

#include "aggregate.h"
#include "combine_axes.h"
#include "extrapolate.h"
#include "filter.h"
#include "helper.h"

#define max(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

dataframe_t mims_unit(dataframe_t *dataframe,
                      int dyanmic_range_low, int dyanmic_range_high,
                      int break_size, time_unit_t time_unit,
                      double noise_level, double k, double spar,
                      double cutoff_low, double cutoff_high,
                      bool allow_truncation);

#endif // _MIMS_UNIT_H_