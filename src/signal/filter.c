#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

static double *cfilter(int nx, double *x, int nf, double *filter, int sides)
{
    int i, j, nshift;
    double z, tmp;
    double *out = malloc(nx * sizeof(double));

    if (sides == 2)
        nshift = nf / 2;
    else
        nshift = 0;

    for (i = 0; i < nx; i++)
    {
        z = 0;
        if (i + nshift - (nf - 1) < 0 || i + nshift >= nx)
        {
            out[i] = NAN;
            continue;
        }
        for (j = MAX(0, nshift + i - nx); j < MIN(nf, i + nshift + 1); j++)
        {
            tmp = x[i + nshift - j];
            if (!isnan(tmp))
                z += filter[j] * tmp;
            else
            {
                out[i] = NAN;
                goto bad;
            }
        }
        out[i] = z;
    bad:
        continue;
    }

    return out;
}

static void rfilter(int nx, double *rx, int nf, double *rf, double *r)
{
    int i, j;
    double sum, tmp;
    // *r = REAL(out), *rx = REAL(x), *rf = REAL(filter);

    for (i = 0; i < nx; i++)
    {
        sum = rx[i];
        for (j = 0; j < nf; j++)
        {
            tmp = r[nf + i - j - 1];
            if (!isnan(tmp))
                sum += tmp * rf[j];
            else
            {
                r[nf + i] = NAN;
                goto bad3;
            }
        }
        r[nf + i] = sum;
    bad3:
        continue;
    }
    return;
}

// method: 0 = convolution, 1 = recursive
static double *stats_filter(int n_x, double *x, int n_filt, double *filt,
                            int method, int sides)
{
    double *y;
    int i = 0;

    if (method == 0) // convolution
    {
        y = cfilter(n_x, x, n_filt, filt, sides);
    }
    else
    {
        double *init = malloc(n_filt * sizeof(double));

        int n_out = n_filt + n_x;
        double *out = malloc(n_out * sizeof(double));
        for (i = 0; i < n_out; i++) // first n_init values are reverse(init)
        {
            if (i < n_filt)
                out[i] = init[n_filt - i - 1];
            else
                out[i] = 0;
        }
        free(init);

        // y = rec_filter(x, filt, out)[-ind]
        rfilter(n_x, x, n_filt, filt, out); // mutates "out"
        int n_y = n_out - n_filt;
        y = malloc(n_y * sizeof(double));
        for (i = 0; i < n_y; i++)
            y[i] = out[i + n_filt];

        free(out);
    }
    return y;
}

double *signal_filter(int n_filt, double *filt, int n_a, double *a, int n_x, double *x)
{
    int i, n_output;
    double *output, *filter_coefs;

    if (n_filt > 0)
    {
        int n_init_x = n_filt - 1;
        int n_concat = n_init_x + n_x;
        double *concat = calloc(n_concat, sizeof(double)); // needs to be calloc
        for (i = n_init_x; i < n_concat; i++)
            concat[i] = x[i - n_init_x];
        filter_coefs = malloc(n_filt * sizeof(double));
        for (i = 0; i < n_filt; i++)
        {
            filter_coefs[i] = filt[i] / a[0];
        }
        double *x1 = stats_filter(n_concat, concat, n_filt, filter_coefs, 0, 1);
        free(filter_coefs);
        free(concat);
        n_output = n_concat;
        for (i = 0; i < n_concat; i++)
            if (isnan(x1[i]))
                n_output--;

        output = malloc(n_output * sizeof(double));
        i = 0;
        int output_i = 0;
        while (output_i <= n_output)
        {
            if (!isnan(x1[i]))
            {
                output[output_i] = x1[i];
                output_i++;
            }
            i++;
        }
        free(x1);
    }

    if (n_a >= 2)
    {
        filter_coefs = malloc((n_a - 1) * sizeof(double));
        for (i = 1; i < n_a; i++)
            filter_coefs[i - 1] = (-1 * a[i]) / a[0];
        double *old_output = output;
        output = stats_filter(n_output, old_output, n_a - 1, filter_coefs, 1, -1);
        free(filter_coefs);
        free(old_output);
    }

    return output;
}
