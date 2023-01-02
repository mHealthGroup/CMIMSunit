#include "butter.h"

transfer_t butter_coeffs(int order, int n_W, double *W, int pass_type, char plane)
{
    // def butter_coeffs(order, W, pass_type="low", plane="z")
    uint8_t stop = 0;
    uint8_t digital = 0;
    int i;
    double T;
    if (pass_type == 1 || pass_type == 2) // stop or high
        stop = 1;
    if (plane == 'z')
        digital = 1;

    // Prewarp to the band edges to s plane
    if (digital)
    {
        T = 2; // sampling frequency of 2 Hz
        // W = (2 / T) * np.tan(np.pi * np.divide(W, T))
        for (i = 0; i < n_W; i++)
        {
            W[i] = (2 / T) * tan(M_PI * (W[i] / T));
        }
    }

    // Generate polelane poles for the prototype butterworth filter
    // source: Kuc
    double C = 1; // default cutoff frequency
    // pole = C * np.exp(1j * np.pi * (2 * np.arange(1, order + 1, 1) + order - 1) / (2 * order))
    // -0.3826834323650897262681+0.9238795325112867384831i
    double complex *pole = malloc(order * sizeof(double complex));
    for (i = 0; i < order; i++)
        pole[i] = C * cexp(I * M_PI * (2 * (i + 1) + order - 1) / (2 * order));

    if (order % 2 == 1)
        pole[(int)((order - 1) / 2)] = -1; // pure real value at exp(i*pi)

    // double *zero;
    zpg_t zpg;
    zpg.zero = calloc(order, sizeof(double));
    zpg.pole = pole;
    zpg.gain = pow(C, order);

    // s-plane frequency transform
    sftrans(order, &zpg, n_W, W, stop);
    order *= 2;
    // Use bilinear transform to convert poles to the z plane
    if (digital)
        bilinear(order, &zpg, T);

    transfer_t tf = as_Arma(order, zpg);
    // a = a.real
    return tf;
}
