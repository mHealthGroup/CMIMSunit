#include <math.h>
#include <stdlib.h>

typedef struct
{
    double ylow;
    double yhigh;
    double f1;
    double f2;
    int kind;
    int na_rm;
} appr_meth;

static double *seq_by_step(double from, double to, double step)
{
    int len = (int)((to - from) / step) + 1;
    double *seq = malloc(len * sizeof(double));
    seq[0] = from;
    seq[len - 1] = to;
    int i;
    for (i = 1; i < len - 1; i++)
    {
        seq[i] = seq[i - 1] + step;
    }

    return seq;
}

static double *seq_by_len(double from, double to, double len)
{
    double step = (to - from) / (len - 1);
    return seq_by_step(from, to, step = step);
}
