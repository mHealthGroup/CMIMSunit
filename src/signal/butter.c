#include "butter.h"

transfer_t butter_coeffs(int order, int n_W, double *W, int pass_type, char plane)
{
    // def butter_coeffs(order, W, pass_type="low", plane="z")
    bool stop = false;
    bool digital = false;
    int i;
    double T;
    if (pass_type == 1 || pass_type == 2) // stop or high
        stop = true;
    if (plane == 'z')
        digital = true;

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
    double complex *pole = calloc(order, sizeof(double complex));
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

// int main(int argc, char **argv)
// {
//     int order = 4;
//     double W[2] = {
//         0.004000000000000000083267,
//         0.100000000000000005551115};
//     int pass_type = 3;

//     transfer_t tf = butter_coeffs(order, 2, W, pass_type, 'z');

//     // $b
//     // 0.0003588407958559478483536  0.0000000000000000000000000 -0.0014353631834237913934144  0.0000000000000000000000000
//     // 0.0021530447751356871985418  0.0000000000000000000000000 -0.0014353631834237913934144  0.0000000000000000000000000
//     // 0.0003588407958559478483536

//     // $a
//     // 1.0000000000000000000000  -7.1989557602839946426343  22.7102400569503402039118 -41.0121387818457563412267
//     // 46.3786038427510902693029 -33.6341446649892930054193  15.2766725132619356486430  -3.9733824787545191092875
//     // 0.4531052730785416482462

//     return 0;
// }
