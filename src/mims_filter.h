#ifndef _FILTER_H_
#define _FILTER_H_

#include <stdlib.h>

#include "mims_unit.h"
#include "signal/filter.h"

dataframe_t *iir(dataframe_t *df, const uint16_t sampling_rate, double *cutoff_freq, const uint8_t order);

#endif // _FILTER_H_