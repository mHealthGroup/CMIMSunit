#ifndef _SIGNAL_BUTTER_H_
#define _SIGNAL_BUTTER_H_

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "arma.h"
#include "bilinear.h"
#include "poly.h"
#include "sftrans.h"

transfer_t butter_coeffs(int order, int n_W, double *W, int pass_type, char plane);

#endif // _SIGNAL_BUTTER_H_