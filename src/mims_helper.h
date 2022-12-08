#ifndef _HELPER_H_
#define _HELPER_H_

#include "mims_unit.h"

#include <math.h>
#include <stdlib.h>

int get_break_size_in_seconds(int break_size, time_unit_t time_unit);
int parse_epoch_string(int break_size, time_unit_t time_unit, int sampling_rate);
int get_sampling_rate(dataframe_t *dataframe);
void segment_data(dataframe_t *dataframe, int break_size, time_unit_t time_unit, double start_time);

int sequence_length(double start, double stop, double step);
double *sequence(double start, double stop, double step);
double *linspace(double start, double stop, int n);

#endif // _HELPER_H_