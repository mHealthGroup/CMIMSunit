#ifndef _HELPER_H_
#define _HELPER_H_

#include "mims_unit.h"

#include <math.h>
#include <stdlib.h>

uint32_t parse_epoch_string(const uint16_t break_size, const time_unit_t time_unit,
                            const uint16_t sampling_rate);
uint16_t get_sampling_rate(const dataframe_t *dataframe);
void segment_data(dataframe_t *dataframe, const uint16_t break_size,
                  const time_unit_t time_unit, const double start_time);

uint32_t sequence_length(const double start, const double stop, const double step);
double *sequence(const double start, const double stop, const double step);
double *linspace(const double start, const double stop, const uint32_t n);

int count_lines(const char *filename);
dataframe_t *read_csv(const char *filename);

dataframe_t *create_dataframe(uint32_t size, double *timestamps, double *x,
                              double *y, double *z, uint32_t n_segments,
                              uint32_t *segments, double *mims_data);
void free_dataframe(dataframe_t *df);
dataframe_t *concat_dataframes(const dataframe_t *df_1, const dataframe_t *df_2);

#endif // _HELPER_H_