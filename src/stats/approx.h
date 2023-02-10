#ifndef _STATS_APPROX_H_
#define _STATS_APPROX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../mims_helper.h"

typedef struct
{
    uint32_t n;
    double *x;
    double *y;
} approx_output_t;

approx_output_t *approx(const uint32_t n, const double *x, const double *y, const uint32_t nout);

#endif // _STATS_APPROX_H_