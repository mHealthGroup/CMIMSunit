#ifndef _AGGREGATE_H_
#define _AGGREGATE_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mims_helper.h"
#include "mims_unit.h"
#include "trapz/trapz.h"

typedef enum angle_unit
{
    degrees,
    radians
} angle_unit_t;

dataframe_t *aggregate(dataframe_t *dataframe, const uint16_t break_size, const time_unit_t time_unit,
                       const uint8_t rectify, const double start_time);

#endif // _AGGREGATE_H_