#ifndef _SFTRANS_H_
#define _SFTRANS_H_

#include <complex.h>
#include <math.h>

#include "arma.h"

void sftrans(int n, zpg_t *zpg, int n_W, double *W, uint8_t stop);

#endif // _SFTRANS_H_