#ifndef _AGGREGATE_H_
#define _AGGREGATE_H_

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "helper.h"
#include "mims_unit.h"
#include "trapz/trapz.h"

typedef enum angle_unit
{
    degrees,
    radians
} angle_unit_t;

dataframe_t aggregate(dataframe_t *dataframe, int break_size, time_unit_t time_unit,
                      bool rectify, double start_time);

#endif // _AGGREGATE_H_