#include "bilinear.h"

void bilinear(int n, zpg_t *zpg, double T)
{
    int i;
    int z = n / 2;
    double *x = calloc(n, sizeof(double));
    for (i = 0; i < n; i++)
        x[i] = (2 - zpg->pole[i] * T) / T;

    double zero_prod;
    for (i = 0; i < z; i++)
    {
        double prod = (2 - zpg->zero[i] * T) / T;
        zero_prod = (i == 0) ? prod : zero_prod * prod;
    }
    double complex pole_prod;
    for (i = 0; i < n; i++)
    {
        double complex prod = (2 - zpg->pole[i] * T) / T;
        pole_prod = (i == 0) ? prod : pole_prod * prod;
    }

    double Zg = creal(zpg->gain * zero_prod / pole_prod);

    double complex *Zp = calloc(n, sizeof(double complex));
    double *Zz = calloc(n, sizeof(double));
    for (i = 0; i < n; i++)
    {
        Zp[i] = (2 + zpg->pole[i] * T) / (2 - zpg->pole[i] * T);
        Zz[i] = (i < z) ? (2 + zpg->zero[i] * T) / (2 - zpg->zero[i] * T) : -1;
    }
    zpg->zero = Zz;
    zpg->pole = Zp;
    zpg->gain = creal(Zg);
}
