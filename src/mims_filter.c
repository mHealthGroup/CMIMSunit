#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mims_filter.h"
#include "mims_helper.h"

// pass_type is always high (2)
// filter type is always butter
dataframe_t iir(dataframe_t *df, int sampling_rate, double *cutoff_freq, int order)
{
  int cutoff_freq_n = 2;
  double nyquist = (double)sampling_rate / 2.0;
  for (int i = 0; i < cutoff_freq_n; i++)
    cutoff_freq[i] = cutoff_freq[i] / nyquist;

  transfer_t coeffs = butter_coeffs(order, cutoff_freq_n, cutoff_freq, 3, 'z');

  dataframe_t result;
  result.size = df->size;
  result.x = signal_filter(coeffs.size, coeffs.b, coeffs.size, coeffs.a, df->size, df->x);
  result.y = signal_filter(coeffs.size, coeffs.b, coeffs.size, coeffs.a, df->size, df->y);
  result.z = signal_filter(coeffs.size, coeffs.b, coeffs.size, coeffs.a, df->size, df->z);

  return result;
}