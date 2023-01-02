#include "poly.h"

double *poly(int n, double *x)
{

    double *y = malloc((n + 1) * sizeof(double));
    double *y_subarr1, *y_subarr2;
    int i, j;

    y[0] = 1;
    for (j = 0; j < n; j++)
    {
        y_subarr1 = malloc((j + 1) * sizeof(double));
        y_subarr2 = malloc((j + 1) * sizeof(double));
        for (i = 1; i <= j + 1; i++)
        {
            y_subarr1[i - 1] = y[i];
            y_subarr2[i - 1] = y[i - 1];
        }
        for (i = 1; i <= j + 1; i++)
            y[i] = y_subarr1[i - 1] - x[j] * y_subarr2[i - 1];
    }

    return y;
}

double complex *complex_poly(int n, double complex *x)
{

    double complex *y = malloc((n + 1) * sizeof(double complex));
    double complex *y_subarr1, *y_subarr2;
    y[0] = 1;
    int i, j;

    y[0] = 1;
    for (j = 0; j < n; j++)
    {
        y_subarr1 = malloc((j + 1) * sizeof(double complex));
        y_subarr2 = malloc((j + 1) * sizeof(double complex));
        for (i = 1; i <= j + 1; i++)
        {
            y_subarr1[i - 1] = y[i];
            y_subarr2[i - 1] = y[i - 1];
        }
        for (i = 1; i <= j + 1; i++)
            y[i] = y_subarr1[i - 1] - x[j] * y_subarr2[i - 1];
    }

    return y;
}
