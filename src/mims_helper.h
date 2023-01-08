#ifndef _HELPER_H_
#define _HELPER_H_

#include "mims_unit.h"

#include <math.h>
#include <stdlib.h>

uint32_t parse_epoch_string(uint16_t break_size, time_unit_t time_unit, uint16_t sampling_rate);
uint16_t get_sampling_rate(dataframe_t *dataframe);
void segment_data(dataframe_t *dataframe, uint16_t break_size, time_unit_t time_unit, double start_time);
uint32_t sequence_length(double start, double stop, double step);
double *sequence(double start, double stop, double step);
double *linspace(double start, double stop, uint32_t n);
dataframe_t concat_dataframes(dataframe_t *df_1, dataframe_t *df_2);
int count_lines(char *filename);
dataframe_t read_csv(char *filename);

#endif // _HELPER_H_