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

static double approx1(const double v, const double *x, const double *y, const uint32_t n,
                      const appr_meth *Meth)
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

static void R_approxfun(const double *x, const double *y, const uint32_t nxy,
                        double *xout, double *yout, const uint32_t nout, const uint8_t method,
                        const double yleft, const double yright, const double f, const uint8_t na_rm)
{
    appr_meth *M = malloc(sizeof(appr_meth)); /* -Wall */
    M->ylow = yleft;
    M->yhigh = yright;
    M->f1 = 1 - f;
    M->f2 = f;
    M->kind = method;
    M->na_rm = na_rm;

    for (uint32_t i = 0; i < nout; i++)
        yout[i] = isnan(xout[i]) ? xout[i] : approx1(xout[i], x, y, nxy, M);

    free(M);
}

static double *C_Approx(const double *x, const double *y, const uint32_t nxy, double *xout, const uint32_t nout, const uint8_t method,
                        const double yleft, const double yright, const double f, const uint8_t na_rm)
{
    double *yout = malloc(nout * sizeof(double));
    R_approxfun(x, y, nxy, xout, yout, nout, method, yleft, yright, f, na_rm);
    return yout;
}

// approx func from
// https://github.com/wch/r-source/blob/79298c499218846d14500255efd622b5021c10ec/src/library/stats/R/approx.R
// "linear" method
approx_output_t *approx(const uint32_t n, const double *x, const double *y, const uint32_t nout)
{
    double *xout = linspace(x[0], x[n - 1], nout);
    double *yout = C_Approx(x, y, n, xout, nout, 1, NAN, NAN, 0, 1);

    approx_output_t *output = malloc(sizeof(approx_output_t));
    output->n = nout;
    output->x = xout;
    output->y = yout;
    return output;
}
