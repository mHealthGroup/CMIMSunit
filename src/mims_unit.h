#ifndef _MIMS_UNIT_H_
#define _MIMS_UNIT_H_

#include <stdlib.h>
#include <stdio.h>

// uint8_t (255)
// uint16_t (65535)
// uint32_t (4294967295)
// uint64_t (18446744073709551615)

typedef struct
{
    uint32_t size;
    double *timestamps; // seconds w/ decimal milliseconds
    double *x;
    double *y;
    double *z;
    uint32_t n_segments;
    uint32_t *segments;
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

dataframe_t *mims_unit_from_filename(char *input_filename,
                                     int8_t dyanmic_range_low, int8_t dyanmic_range_high,
                                     uint16_t break_size, time_unit_t time_unit,
                                     float noise_level, float k, float spar,
                                     float cutoff_low, float cutoff_high,
                                     uint8_t allow_truncation);
dataframe_t *mims_unit(dataframe_t *dataframe,
                       int8_t dyanmic_range_low, int8_t dyanmic_range_high,
                       uint16_t break_size, time_unit_t time_unit,
                       float noise_level, float k, float spar,
                       float cutoff_low, float cutoff_high,
                       uint8_t allow_truncation);
dataframe_t *custom_mims_unit(dataframe_t *dataframe,
                              int8_t dyanmic_range_low, int8_t dyanmic_range_high,
                              uint16_t break_size, time_unit_t time_unit,
                              float noise_level, float k, float spar,
                              float cutoff_low, float cutoff_high,
                              uint8_t allow_truncation);
dataframe_t *custom_mims_unit_before_after_dataframe(dataframe_t *dataframe,
                                                     int8_t dyanmic_range_low, int8_t dyanmic_range_high,
                                                     uint16_t break_size, time_unit_t time_unit,
                                                     float noise_level, float k, float spar,
                                                     float cutoff_low, float cutoff_high,
                                                     uint8_t allow_truncation,
                                                     dataframe_t *before_df, dataframe_t *after_df);

#endif // _MIMS_UNIT_H_
