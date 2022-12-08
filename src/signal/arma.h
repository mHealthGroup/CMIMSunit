#ifndef _ARMA_H_
#define _ARMA_H_

#include <complex.h>
#include <stdlib.h>

#include "poly.h"

typedef struct
{
    int size;
    double *a;
    double *b;
} transfer_t;

typedef struct
{
    double *zero;
    double complex *pole;
    double gain;
} zpg_t;

transfer_t as_Arma(int n, zpg_t zpg);

#endif // _ARMA_H_