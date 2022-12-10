#include <stdio.h>
#include <stdlib.h>
#include "spline.h"

static void natural_spline(int n, double *x, double *y, double *b, double *c, double *d)
{
    if (n < 2)
    {
        printf("natural_spline n < 2");
        return;
    }

    x--;
    y--;
    b--;
    c--;
    d--;

    if (n < 3)
    {
        double t = (y[2] - y[1]);
        b[1] = t / (x[2] - x[1]);
        b[2] = b[1];
        c[1] = c[2] = d[1] = d[2] = 0.0;
        return;
    }

    const int nm1 = n - 1;
    int i;

    /* Set up the tridiagonal system */
    /* b = diagonal, d = offdiagonal, c = right hand side */

    d[1] = x[2] - x[1];
    c[2] = (y[2] - y[1]) / d[1];
    for (i = 2; i < n; i++)
    {
        d[i] = x[i + 1] - x[i];
        b[i] = 2.0 * (d[i - 1] + d[i]);
        c[i + 1] = (y[i + 1] - y[i]) / d[i];
        c[i] = c[i + 1] - c[i];
    }

    /* Gaussian elimination */

    for (i = 3; i < n; i++)
    {
        double t = d[i - 1] / b[i - 1];
        b[i] = b[i] - t * d[i - 1];
        c[i] = c[i] - t * c[i - 1];
    }

    /* Backward substitution */

    c[nm1] = c[nm1] / b[nm1];
    for (i = n - 2; i > 1; i--)
        c[i] = (c[i] - d[i] * c[i + 1]) / b[i];

    /* End conditions */

    c[1] = c[n] = 0.0;

    /* Get cubic coefficients */

    b[1] = (y[2] - y[1]) / d[1] - d[i] * c[2];
    c[1] = 0.0;
    d[1] = c[2] / d[1];
    b[n] = (y[n] - y[nm1]) / d[nm1] + d[nm1] * c[nm1];
    for (i = 2; i < n; i++)
    {
        b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i]);
        d[i] = (c[i + 1] - c[i]) / d[i];
        c[i] = 3.0 * c[i];
    }
    c[n] = 0.0;
    d[n] = 0.0;

    return;
}

static void fmm_spline(int n, double *x, double *y, double *b, double *c, double *d)
{
    /* Adjustment for 1-based arrays */
    x--;
    y--;
    b--;
    c--;
    d--;

    if (n < 2)
    {
        return;
    }

    if (n < 3)
    {
        double t = (y[2] - y[1]);
        b[1] = t / (x[2] - x[1]);
        b[2] = b[1];
        c[1] = c[2] = d[1] = d[2] = 0.0;
        return;
    }

    const int nm1 = n - 1;
    int i;

    /* Set up tridiagonal system */
    /* b = diagonal, d = offdiagonal, c = right hand side */

    d[1] = x[2] - x[1];
    c[2] = (y[2] - y[1]) / d[1]; /* = +/- Inf	for x[1]=x[2] -- problem? */
    for (i = 2; i < n; i++)
    {
        d[i] = x[i + 1] - x[i];
        b[i] = 2.0 * (d[i - 1] + d[i]);
        c[i + 1] = (y[i + 1] - y[i]) / d[i];
        c[i] = c[i + 1] - c[i];
    }

    /* End conditions. */
    /* Third derivatives at x[0] and x[n-1] obtained */
    /* from divided differences */

    b[1] = -d[1];
    b[n] = -d[nm1];
    c[1] = c[n] = 0.0;
    if (n > 3)
    {
        c[1] = c[3] / (x[4] - x[2]) - c[2] / (x[3] - x[1]);
        c[n] = c[nm1] / (x[n] - x[n - 2]) - c[n - 2] / (x[nm1] - x[n - 3]);
        c[1] = c[1] * d[1] * d[1] / (x[4] - x[1]);
        c[n] = -c[n] * d[nm1] * d[nm1] / (x[n] - x[n - 3]);
    }

    /* Gaussian elimination */

    for (i = 2; i <= n; i++)
    {
        double t = d[i - 1] / b[i - 1];
        b[i] = b[i] - t * d[i - 1];
        c[i] = c[i] - t * c[i - 1];
    }

    /* Backward substitution */

    c[n] = c[n] / b[n];
    for (i = nm1; i >= 1; i--)
        c[i] = (c[i] - d[i] * c[i + 1]) / b[i];

    /* c[i] is now the sigma[i-1] of the text */
    /* Compute polynomial coefficients */

    b[n] = (y[n] - y[n - 1]) / d[n - 1] + d[n - 1] * (c[n - 1] + 2.0 * c[n]);
    for (i = 1; i <= nm1; i++)
    {
        b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i]);
        d[i] = (c[i + 1] - c[i]) / d[i];
        c[i] = 3.0 * c[i];
    }
    c[n] = 3.0 * c[n];
    d[n] = d[nm1];
    return;
}

static void spline_coef(int method, int n, double *x, double *y,
                        double *b, double *c, double *d)
{
    switch (method)
    {
    case 1:
        // periodic_spline(n, x, y, b, c, d);
        // never used?
        break;

    case 2:
        natural_spline(n, x, y, b, c, d);
        break;

    case 3:
        fmm_spline(n, x, y, b, c, d);
        break;
    }
}

// TODO Create struct for return type (equivalent of python dict)
static Z_struct_t SplineCoef(int method, int n, double *x, int m, double *y)
{
    // x = PROTECT(coerceVector(x, REALSXP));
    // y = PROTECT(coerceVector(y, REALSXP));
    // R_xlen_t n = XLENGTH(x);
    // int method = asInteger(method);

    if (m != n)
    {
        printf("Error: inputs of different lengths");
    }

    // ans, nm;
    // b = PROTECT(allocVector(REALSXP, n));
    // c = PROTECT(allocVector(REALSXP, n));
    // d = PROTECT(allocVector(REALSXP, n));
    // double *rb = REAL(b), *rc = REAL(c), *rd = REAL(d);

    double *b = malloc(n * sizeof(double));
    double *c = malloc(n * sizeof(double));
    double *d = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
    {
        b[i] = c[i] = d[i] = 0;
    }

    spline_coef(method, n, x, y, b, c, d);

    FILE *fpt;
    fpt = fopen("/Users/arytonhoi/Kode/mhealth/pymims/planning/C/spline/spline_coef.csv", "w+");
    fprintf(fpt, "b, c, d\n");
    for (int i = 0; i < n; i++)
    {
        fprintf(fpt, "%.10f, %.10f, %.10f\n", b[i], c[i], d[i]);
    }
    fclose(fpt);

    // // ans = PROTECT(allocVector(VECSXP, 7));
    // // SET_VECTOR_ELT(ans, 0, ScalarInteger(method));
    // // SET_VECTOR_ELT(ans, 1, (n > INT_MAX) ? ScalarReal((double)n) : ScalarInteger((int)n));
    // // SET_VECTOR_ELT(ans, 2, x);
    // // SET_VECTOR_ELT(ans, 3, y);
    // // SET_VECTOR_ELT(ans, 4, b);
    // // SET_VECTOR_ELT(ans, 5, c);
    // // SET_VECTOR_ELT(ans, 6, d);
    Z_struct_t ans;
    ans.method = method;
    ans.n = n;
    ans.x = x;
    ans.y = y;
    ans.b = b;
    ans.c = c;
    ans.d = d;

    // // nm = allocVector(STRSXP, 7);
    // // setAttrib(ans, R_NamesSymbol, nm);
    // // SET_STRING_ELT(nm, 0, mkChar("method"));
    // // SET_STRING_ELT(nm, 1, mkChar("n"));
    // // SET_STRING_ELT(nm, 2, mkChar("x"));
    // // SET_STRING_ELT(nm, 3, mkChar("y"));
    // // SET_STRING_ELT(nm, 4, mkChar("b"));
    // // SET_STRING_ELT(nm, 5, mkChar("c"));
    // // SET_STRING_ELT(nm, 6, mkChar("d"));
    // // UNPROTECT(6);

    return ans;
}

static void spline_eval(int method, int nu, double *u, double *v,
                        int n, double *x, double *y, double *b, double *c, double *d)
{
    /* Evaluate  v[l] := spline(u[l], ...),	    l = 1,..,nu, i.e. 0:(nu-1)
     * Nodes x[i], coef (y[i]; b[i],c[i],d[i]); i = 1,..,n , i.e. 0:(*n-1)
     */
    const int n_1 = n - 1;
    int i, l;
    double dx;

    // if (method == 1 && n > 1)
    // { /* periodic */
    //     dx = x[n_1] - x[0];
    //     for (l = 0; l < nu; l++)
    //     {
    //         v[l] = fmod(u[l] - x[0], dx);
    //         if (v[l] < 0.0)
    //             v[l] += dx;
    //         v[l] += x[0];
    //     }
    // }
    // else

    // commenting out because errors
    // for (l = 0; l < nu; l++)
    //     v[l] = u[l];

    i = 0;
    for (l = 0; l < nu; l++)
    {
        double ul = u[l];
        if (ul < x[i] || (i < n_1 && x[i + 1] < ul))
        {
            /* reset i  such that  x[i] <= ul <= x[i+1] : */
            i = 0;
            int j = n;
            do
            {
                int k = (i + j) / 2;
                if (ul < x[k])
                    j = k;
                else
                    i = k;
            } while (j > i + 1);
        }
        dx = ul - x[i];
        /* for natural splines extrapolate linearly left */
        double tmp = (method == 2 && ul < x[0]) ? 0.0 : d[i];

        v[l] = y[i] + dx * (b[i] + dx * (c[i] + dx * tmp));
    }
}

static double *SplineEval(int nu, double *xout, Z_struct_t z)
{
    // xout = PROTECT(coerceVector(xout, REALSXP));
    // R_xlen_t nu = XLENGTH(xout)

    // nx = asXlen(getListElement(z, "n"));
    int nx = z.n;

    // SEXP yout = PROTECT(allocVector(REALSXP, nu));
    double *yout = malloc(nu * sizeof(double));

    // int method = asInteger(getListElement(z, "method"));
    // SEXP x = getListElement(z, "x"), y = getListElement(z, "y"),
    //      b = getListElement(z, "b"), c = getListElement(z, "c"),
    //      d = getListElement(z, "d");
    // spline_eval(method, nu, REAL(xout), REAL(yout),
    //             nx, REAL(x), REAL(y), REAL(b), REAL(c), REAL(d));
    spline_eval(z.method, nu, xout, yout, nx, z.x, z.y, z.b, z.c, z.d);

    // UNPROTECT(2);
    return yout;
}

// Method 1: periodic (unimplemented because unused), 2: natural, 3: fmm
double *Spline(int x_len, double *x, int y_len, double *y,
               int xout_len, double *xout, int method)
{
    Z_struct_t z = SplineCoef(method, x_len, x, y_len, y);
    double *yout = SplineEval(xout_len, xout, z);

    // FILE *fpt;
    // fpt = fopen("/Users/arytonhoi/Kode/mhealth/pymims/planning/C/spline/spline_eval.csv", "w+");
    // fprintf(fpt, "x, y\n");
    // for (int i = 0; i < xout_len; i++)
    // {
    //     fprintf(fpt, "%f, %.15f\n", xout[i], yout[i]);
    // }
    // fclose(fpt);

    return yout;
}

// int main(int argc, char **argv)
// {
//     return 0;
// }
