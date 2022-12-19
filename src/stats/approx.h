#ifndef _STATS_APPROX_H_
#define _STATS_APPROX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../mims_helper.h"

typedef struct
{
    int n;
    double *x;
    double *y;
} approx_output_t;

approx_output_t approx(uint32_t n, double *x, double *y, uint32_t nout);

#endif // _STATS_APPROX_H_