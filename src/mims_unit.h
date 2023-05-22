#ifndef _MIMS_UNIT_H_
#define _MIMS_UNIT_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

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

dataframe_t *mims_unit_from_filename(const char *input_filename,
                                     const int8_t dyanmic_range_low, const int8_t dyanmic_range_high,
                                     const uint16_t break_size, const time_unit_t time_unit,
                                     const float noise_level, const float k, const float spar,
                                     const float cutoff_low, const float cutoff_high,
                                     const uint8_t allow_truncation);
dataframe_t *mims_unit(dataframe_t *dataframe,
                       const int8_t dyanmic_range_low, const int8_t dyanmic_range_high,
                       const uint16_t break_size, const time_unit_t time_unit,
                       const float noise_level, const float k, const float spar,
                       const float cutoff_low, const float cutoff_high,
                       const uint8_t allow_truncation);
dataframe_t *custom_mims_unit_before_after_dataframe(dataframe_t *dataframe,
                                                     const int8_t dyanmic_range_low, const int8_t dyanmic_range_high,
                                                     const uint16_t break_size, const time_unit_t time_unit,
                                                     const float noise_level, const float k, const float spar,
                                                     const float cutoff_low, const float cutoff_high,
                                                     const uint8_t allow_truncation,
                                                     const dataframe_t *before_df, const dataframe_t *after_df);
dataframe_t *custom_mims_unit(const dataframe_t *dataframe,
                              const int8_t dyanmic_range_low, const int8_t dyanmic_range_high,
                              const uint16_t break_size, const time_unit_t time_unit,
                              const float noise_level, const float k, const float spar,
                              const float cutoff_low, const float cutoff_high,
                              const uint8_t allow_truncation);

#endif // _MIMS_UNIT_H_
