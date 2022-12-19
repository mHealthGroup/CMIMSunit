#include "trapz.h"

double catools_trapz(uint16_t n, double *x, double *y)
{
    // computes the integral of y with respect to x using trapezoidal integration.

    double output = 0.0;
    for (uint16_t i = 1; i < n; i++)
        output += (x[i] - x[i - 1]) * (y[i] + y[i - 1]);

    return output / 2.0;
}
