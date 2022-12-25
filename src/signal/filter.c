#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

double *cfilter(int nx, double *x, int nf, double *filter, int sides)
{
    int i, j, nshift;
    double z, tmp;
    double *out = calloc(nx, sizeof(double));

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

void rfilter(int nx, double *rx, int nf, double *rf, double *r)
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
double *stats_filter(int n_x, double *x, int n_filt, double *filt,
                     int method, int sides)
{
    double *y;
    int i = 0;
    //     method <- match.arg(method)
    //     x <- as.ts(x)
    //     storage.mode(x) <- "double"
    //     xtsp <- tsp(x)
    //     n <- as.integer(NROW(x))
    //     if (is.na(n)) stop(gettextf("invalid value of %s", "NROW(x)"), domain = NA)
    //     nser <- NCOL(x)
    //     filter <- as.double(filter)
    //     nfilt <- as.integer(length(filter))
    //     if (is.na(nfilt)) stop(gettextf("invalid value of %s", "length(filter)"),
    //                            domain = NA)
    //     if(anyNA(filter)) stop("missing values in 'filter'")

    if (method == 0) // convolution
    {
        y = cfilter(n_x, x, n_filt, filt, sides);
    }
    else
    {
        double *init = calloc(n_filt, sizeof(double));
        // int *ind = calloc(n_filt, sizeof(int));
        // for (i = 0; i < n_filt; i++)
        //     ind[i] = i;

        int n_out = n_filt + n_x;
        double *out = calloc(n_out, sizeof(double));
        for (i = 0; i < n_out; i++) // first n_init values are reverse(init)
        {
            if (i < n_filt)
                out[i] = init[n_filt - i - 1];
            else
                out[i] = 0;
        }

        // y = rec_filter(x, filt, out)[-ind]
        rfilter(n_x, x, n_filt, filt, out); // mutates "out"
        int n_y = n_out - n_filt;
        y = calloc(n_y, sizeof(double));
        for (i = 0; i < n_y; i++)
            y[i] = out[i + n_filt];
    }
    return y;
}

double *signal_filter(int n_filt, double *filt, int n_a, double *a, int n_x, double *x)
{
    // int n_init = n_filt - 1;
    int i, n_output;
    // double *init = calloc(n_init, sizeof(double));
    double *output, *filter_coefs;

    if (n_filt > 0)
    {
        // double *x1 = stats_filter(
        //     x = np.concatenate((init_x, x)),
        //     filter_coeffs = filt / a[0],
        //     method = "convolution",
        //     sides = 1, )
        // output = x1[~np.isnan(x1)]

        // x1 <- stats::filter(c(init.x, x), filt / a[1], sides = 1)
        // if(all(is.na(x1)))
        //     return(x)
        // x <- na.omit(x1, filt / a[1] , sides = 1)

        int n_init_x = n_filt - 1;
        int n_concat = n_init_x + n_x;
        double *concat = calloc(n_concat, sizeof(double));
        for (i = n_init_x; i < n_concat; i++)
            concat[i] = x[i - n_init_x];
        filter_coefs = calloc(n_filt, sizeof(double));
        for (i = 0; i < n_filt; i++)
        {
            filter_coefs[i] = filt[i] / a[0];
        }
        double *x1 = stats_filter(n_concat, concat, n_filt, filter_coefs, 0, 1);
        n_output = n_concat;
        for (i = 0; i < n_concat; i++)
            if (isnan(x1[i]))
                n_output--;

        output = calloc(n_output, sizeof(double));
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
    }

    if (n_a >= 2)
    {
        // output = stats_filter(
        //     x=x, filter_coeffs=np.divide(-a[1:], a[0]), method="recursive", init=init);

        // x <- stats::filter(x, -a[-1] / a[1], method = "recursive", init = init)
        filter_coefs = calloc(n_a - 1, sizeof(double));
        for (i = 1; i < n_a; i++)
            filter_coefs[i - 1] = (-1 * a[i]) / a[0];
        output = stats_filter(n_output, output, n_a - 1, filter_coefs, 1, -1);
    }

    return output;
}

// int main(int argc, char **argv)
// {
//     double filt[9] = {
//         0.0003588407958559478483536,
//         0.0000000000000000000000000,
//         -0.0014353631834237913934144,
//         0.0000000000000000000000000,
//         0.0021530447751356871985418,
//         0.0000000000000000000000000,
//         -0.0014353631834237913934144,
//         0.0000000000000000000000000,
//         0.0003588407958559478483536,

//     };

//     double a[9] = {
//         1.0000000000000000000000,
//         -7.1989557602839946426343,
//         22.7102400569503402039118,
//         -41.0121387818457563412267,
//         46.3786038427510902693029,
//         -33.6341446649892930054193,
//         15.2766725132619356486430,
//         -3.9733824787545191092875,
//         0.4531052730785416482462,
//     };

//     double x[21] = {
//         -0.2540000000000000035527,
//         -0.2431731268070462803621,
//         -0.2242802170493227720272,
//         -0.2345699232520407906399,
//         -0.2410431758077686004160,
//         -0.1990000000000000102141,
//         -0.2087867856385630105365,
//         -0.2200292455890628939841,
//         -0.2320821309642192020739,
//         -0.2594364509687815401051,
//         -0.2770000000000000239808,
//         -0.2688028142034358247692,
//         -0.2417369317174523912772,
//         -0.2154482763349978013956,
//         -0.2136492266466666067881,
//         -0.2340000000000000135447,
//         -0.2329666319280039032957,
//         -0.2243971288436006350508,
//         -0.2254370715060204088953,
//         -0.2556853034408685387824,
//         -0.3009999999999999897859,
//     };

//     // [1] -9.114556214741075619232e-05 -7.434133079996452708391e-04
//     // [3] -2.997760562393522298236e-03 -8.170856931057240632454e-03
//     // [5] -1.731498733777252913013e-02 -3.087768831524188850590e-02
//     // [7] -4.859087168207139317833e-02 -6.951725153958990266467e-02
//     // [9] -9.226469706647616453310e-02 -1.152894869852837400614e-01
//     // [11] -1.371720607051950258093e-01 -1.568003872775814711016e-01
//     // [13] -1.734146886905368889487e-01 -1.865211877674857743337e-01
//     // [15] -1.957746011565263677401e-01 -2.009455092585747393308e-01
//     // [17] -2.019922889843708269098e-01 -1.991367680539182782873e-01
//     // [19] -1.928505058085419332503e-01 -1.837867141367923173867e-01
//     // [21] -1.727524251337742844381e-01

//     double *output = signal_filter(9, filt, 9, a, 21, x);
//     return 0;
// }
