#include "arma.h"

transfer_t as_Arma(int n, zpg_t zpg)
{
    transfer_t tf;
    tf.size = n + 1;
    tf.a = calloc(n, sizeof(double));
    tf.b = calloc(n, sizeof(double));
    int i;
    tf.b = poly(n, zpg.zero);
    for (i = 0; i < tf.size; i++)
        tf.b[i] = creal(tf.b[i] * zpg.gain);

    double complex *complex_a = complex_poly(n, zpg.pole);
    for (i = 0; i < tf.size; i++)
        tf.a[i] = creal(complex_a[i]);

    return tf;
}
