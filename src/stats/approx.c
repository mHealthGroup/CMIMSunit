#include "stats_approx.h"

typedef struct
{
    double ylow;
    double yhigh;
    double f1;
    double f2;
    int kind;
    int na_rm;
} appr_meth;

static double approx1(double v, double *x, double *y, int n,
                      appr_meth *Meth)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */

    if (!n)
        return NAN;

    int i = 0;
    int j = n - 1;
    /* handle out-of-domain points */
    if (v < x[i])
        return Meth->ylow;
    if (v > x[j])
        return Meth->yhigh;

    /* find the correct interval by bisection */
    while (i < j - 1)
    { /* x[i] <= v <= x[j] */
        int ij = (i + j) / 2;
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

static void R_approxfun(double *x, double *y, int nxy,
                        double *xout, double *yout, int nout, int method,
                        double yleft, double yright, double f, int na_rm)
{
    appr_meth M = {0.0, 0.0, 0.0, 0.0, 0}; /* -Wall */

    M.f2 = f;
    M.f1 = 1 - f;
    M.kind = method;
    M.ylow = yleft;
    M.yhigh = yright;
    M.na_rm = na_rm;
    int i;
    for (i = 0; i < nout; i++)
        yout[i] = isnan(xout[i]) ? xout[i] : approx1(xout[i], x, y, nxy, &M);
}

static double *C_Approx(double *x, double *y, int nxy, double *xout, int nout, int method,
                        double yleft, double yright, double f, int na_rm)
{
    double *yout = calloc(nout, sizeof(double));
    R_approxfun(x, y, nxy, xout, yout, nout, method, yleft, yright, f, na_rm);
    return yout;
}

// approx func from
// https://github.com/wch/r-source/blob/79298c499218846d14500255efd622b5021c10ec/src/library/stats/R/approx.R
// "linear" method
approx_output_t approx(int n, double *x, double *y, int nout)
{
    double *xout = linspace(x[0], x[n - 1], nout);
    double *yout = C_Approx(x, y, n, xout, nout, 1, NAN, NAN, 0, 1);

    approx_output_t output = {.x = xout, .y = yout};
    return output;
}

// int main(int argc, char **argv)
// {
//     // double sub_t[5] = {
//     //     1495984789.299999952316,
//     //     1495984789.30999994278,
//     //     1495984789.319999933243,
//     //     1495984789.329999923706,
//     //     1495984789.339999914169,
//     // };
//     // double sub_value[5] = {
//     //     1.870999999999999996447,
//     //     3.470134220599003782581,
//     //     5.674122271358177371781,
//     //     7.305276141552649704636,
//     //     7.917618306893985824502,
//     // };
//     // $x
//     // [1] 1495984789.299999952316 1495984789.309999942780{ 1495984789.319999933243
//     // [4] 1495984789.329999923706 1495984789.339999914169

//     // $y
//     // [1] 1.870999999999999996447 3.470134220599003782581 5.674122271358177371781
//     // [4] 7.305276141552649704636 7.917618306893985824502}

//     double sub_t[5] = {
//         1495984789.349999904633,
//         1495984789.359999895096,
//         1495984789.369999885559,
//         1495984789.380000114441,
//         1495984789.390000104904,
//     };
//     double sub_value[5] = {
//         7.996000000000000440536,
//         7.099310714444851733163,
//         5.386598321637007913409,
//         3.708755456133160155474,
//         2.520687807014621473201,
//     };

//     // $x
//     // [1] 1495984789.349999904633 1495984789.359999895096 1495984789.369999885559
//     // [4] 1495984789.380000114441 1495984789.390000104904

//     // $y
//     // [1] 7.996000000000000440536 7.099310714444851733163 5.386598321637007913409
//     // [4] 3.708755456133160155474 2.520687807014621473201

//     int n_over = 5;
//     int method = 1;
//     double yleft = 0;
//     double yright = 0;
//     double f = 0;
//     int na_rm = 1;

//     // approx_output_t thing = approx(5, sub_t, sub_value, n_over);
//     approx_output_t thing = approx(5, sub_t, sub_value, 20);
//     return 0;
// }