#include "sftrans.h"

void sftrans(int n, zpg_t *zpg, int n_W, double *W, bool stop)
{
    int C = 1;
    int i;
    double complex *b = calloc(n, sizeof(double complex));
    double complex temp;
    int p = n;
    int z = 0;

    if (n_W == 2)
    {
        double Fl = W[0];
        double Fh = W[1];
        if (stop)
        {
            // Not debugged
            double zero_prod = 1;
            double pole_prod = 1;
            for (i = 0; i < n; i++)
            {
                zero_prod = zero_prod * (-1 * zpg->zero[i]);
                pole_prod = pole_prod * (-1 * zpg->pole[i]);
            }

            zpg->gain = zpg->gain * (zero_prod / pole_prod);
            double b_numerator = C * (Fh - Fl) / 2;
            for (i = 0; i < n; i++)
                b[i] = b_numerator / zpg->pole[i];
            zpg->pole = calloc(2 * n, sizeof(double complex));
            for (i = 0; i < 2 * n; i++)
            {
                temp = sqrt(pow(b[i], 2) - Fh * Fl) * I;
                zpg->pole[i] = (i < n) ? b[i % n] + temp : b[i % n] - temp;
            }

            // double complex extend_value = sqrt(I * -1 * Fh * Fl);
            // double complex extend[2] = {extend_value, -1 * extend_value};

            b = calloc(n, sizeof(double complex));
            for (i = 0; i < n; i++)
                b[i] = (C * (Fh - Fl) / 2) / zpg->zero[i];

            zpg->zero = calloc(2 * n, sizeof(double));
            for (i = 0; i < 2 * n; i++)
            {
                temp = sqrt(I * pow(b[i % n], 2) - Fh * Fl);
                // probably have complex type mismatch since b is complex and zero isn't
                zpg->zero[i] = (i < n) ? b[i % n] + temp : b[i % n] - temp;
            }
        }
        else
        {
            // gain = gain * (C / (Fh - Fl)) ** (z - p) // z == p
            zpg->gain = zpg->gain * pow(C / (Fh - Fl), z - p);
            for (i = 0; i < n; i++)
                b[i] = zpg->pole[i] * (Fh - Fl) / (2 * C);

            zpg->pole = calloc(2 * n, sizeof(double complex));
            for (i = 0; i < 2 * n; i++)
            {
                temp = csqrt(cpow(b[i % n], 2) - Fh * Fl);
                zpg->pole[i] = (i < n) ? b[i % n] + temp : b[i % n] - temp;
            }

            // for (i = 0; i < n; i++)
            //     b[i] = zpg->zero[i] * (Fh - Fl) / (2 * C);
            // zpg->zero = calloc(2 * n, sizeof(double complex));
            // for (i = 0; i < 2 * n; i++)
            // {
            //     temp = sqrt(I * pow(b[i % n], 2) - Fh * Fl);
            //     zpg->zero[i] = (i < n) ? b[i % n] + temp : b[i % n] - temp;
            // }
            zpg->zero = calloc(p, sizeof(double));
        }
    }
    else
    {
        // not debugged
        double *Fc = W;
        if (stop)
        {
            double zero_prod = 1;
            double pole_prod = 1;
            for (i = 0; i < n; i++)
            {
                zero_prod = zero_prod * (-1 * zpg->zero[i]);
                pole_prod = pole_prod * (-1 * zpg->pole[i]);
            }
            zpg->gain = zpg->gain * (zero_prod / pole_prod);
            for (i = 0; i < n; i++)
                zpg->zero[i] = C * Fc[i] / zpg->zero[i];
        }
        else
        {
            // gain = gain * (C / Fc) ** (z - p) // z == p
            for (i = 0; i < n; i++)
            {
                zpg->pole[i] = Fc[i] * zpg->pole[i] / C;
                zpg->zero[i] = Fc[i] * zpg->zero[i] / C;
            }
        }
    }
}
