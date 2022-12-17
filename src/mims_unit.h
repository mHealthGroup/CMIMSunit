#ifndef _MIMS_UNIT_H_
#define _MIMS_UNIT_H_

#include <stdlib.h>
#include <stdio.h>

typedef struct
{
    int size;
    double *timestamps; // seconds w/ decimal milliseconds
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

#include "mims_aggregate.h"
#include "mims_combine_axes.h"
#include "mims_extrapolate.h"
#include "mims_filter.h"
#include "mims_helper.h"

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
                      uint8_t allow_truncation);

#endif // _MIMS_UNIT_H_