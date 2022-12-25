#include "approx.h"

typedef struct
{
    double ylow;
    double yhigh;
    double f1;
    double f2;
    uint8_t kind;
    uint8_t na_rm;
} appr_meth;

static double approx1(double v, double *x, double *y, uint32_t n,
                      appr_meth *Meth)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */

    if (!n)
        return NAN;

    uint32_t i = 0;
    uint32_t j = n - 1;
    /* handle out-of-domain points */
    if (v < x[i])
        return Meth->ylow;
    if (v > x[j])
        return Meth->yhigh;

    /* find the correct interval by bisection */
    while (i < j - 1)
    { /* x[i] <= v <= x[j] */
        uint32_t ij = (i + j) / 2;
        /* i+1 <= ij <= j-1 */
        if (v < x[ij])
            j = ij;
        else
            i = ij;
        /* still i < j */
    }
    /* provably have i == j-1 */

    /* interpolation */

    if (v == x[j])
        return y[j];
    if (v == x[i])
        return y[i];
    /* impossible: if(x[j] == x[i]) return y[i]; */

    if (Meth->kind == 1) /* linear */
        return y[i] + (y[j] - y[i]) * ((v - x[i]) / (x[j] - x[i]));
    else /* 2 : constant */
        return (Meth->f1 != 0.0 ? y[i] * Meth->f1 : 0.0) + (Meth->f2 != 0.0 ? y[j] * Meth->f2 : 0.0);
} /* approx1() */

static void R_approxfun(double *x, double *y, uint32_t nxy,
                        double *xout, double *yout, uint32_t nout, uint8_t method,
                        double yleft, double yright, double f, uint8_t na_rm)
{
    appr_meth M = {0.0, 0.0, 0.0, 0.0, 0}; /* -Wall */

    M.f2 = f;
    M.f1 = 1 - f;
    M.kind = method;
    M.ylow = yleft;
    M.yhigh = yright;
    M.na_rm = na_rm;
    uint32_t i;
    for (i = 0; i < nout; i++)
        yout[i] = isnan(xout[i]) ? xout[i] : approx1(xout[i], x, y, nxy, &M);
}

static double *C_Approx(double *x, double *y, uint32_t nxy, double *xout, uint32_t nout, uint8_t method,
                        double yleft, double yright, double f, uint8_t na_rm)
{
    double *yout = malloc(nout * sizeof(double));
    R_approxfun(x, y, nxy, xout, yout, nout, method, yleft, yright, f, na_rm);
    return yout;
}

// approx func from
// https://github.com/wch/r-source/blob/79298c499218846d14500255efd622b5021c10ec/src/library/stats/R/approx.R
// "linear" method
approx_output_t approx(uint32_t n, double *x, double *y, uint32_t nout)
{
    double *xout = linspace(x[0], x[n - 1], nout);
    double *yout = C_Approx(x, y, n, xout, nout, 1, NAN, NAN, 0, 1);

    approx_output_t output = {.n = nout, .x = xout, .y = yout};
    return output;
}
