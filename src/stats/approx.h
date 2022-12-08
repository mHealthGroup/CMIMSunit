#ifndef _STATS_APPROX_H_
#define _STATS_APPROX_H_

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "../helper.h"

typedef struct
{
    int n;
    double *x;
    double *y;
} approx_output_t;

approx_output_t approx(int n, double *x, double *y, int nout);

#endif // _STATS_APPROX_H_