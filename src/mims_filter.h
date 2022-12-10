#ifndef _FILTER_H_
#define _FILTER_H_

#include <stdlib.h>

#include "mims_unit.h"
#include "signal/butter.h"
#include "signal/filter.h"

dataframe_t iir(dataframe_t *df, int sampling_rate, double *cutoff_freq, int order);

#endif // _FILTER_H_