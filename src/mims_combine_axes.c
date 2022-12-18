#include "mims_combine_axes.h"

double *sum_up(dataframe_t *dataframe)
{
  uint32_t n = dataframe->size;
  double *output = malloc(n * sizeof(double));
  for (uint32_t i = 0; i < n; i++)
    output[i] = dataframe->x[i] + dataframe->y[i] + dataframe->z[i];

  return output;
}
