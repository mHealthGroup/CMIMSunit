#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mims_filter.h"
#include "mims_helper.h"

typedef struct
{
  int size;
  double *a;
  double *b;
} transfer_t;

// pass_type is always high (2)
// filter type is always butter
dataframe_t iir(dataframe_t *df, uint16_t sampling_rate, double *cutoff_freq, uint8_t order)
{
  uint8_t cutoff_freq_n = 2;
  float nyquist = (float)sampling_rate / 2.0;
  for (uint8_t i = 0; i < cutoff_freq_n; i++)
    cutoff_freq[i] = cutoff_freq[i] / nyquist;

  // transfer_t coeffs = butter_coeffs(order, cutoff_freq_n, cutoff_freq, 3, 'z');
  //
  // butter_coeffs function requires complex.h, which is problematic for some compilers. To
  // temporarily workaround this issue, we're trading flexibility for simplicity by hardcoding
  // butter_coeffs inputs as:
  //
  //    order = 4
  //    cutoff_freq_n = 2
  //    cuttoff_freq = [
  //      0.0062832679918897721,
  //      0.15838444032453627
  //    ]
  //
  // Note sampling rate always = 100hz after our resampling step, so cutoff_freq doesn't change.
  // With these hardcoded inputs, the pre-computed butter coefficients are:

  transfer_t coeffs = {.size = 9};
  double coeffs_a[9] =
      {1.0,
       -7.1989557602839946,
       22.71024005695034,
       -41.012138781845756,
       46.37860384275109,
       -33.634144664989293,
       15.276672513261936,
       -3.9733824787545187,
       0.4531052730785417};
  double coeffs_b[9] =
      {0.00035884079585594785,
       0,
       -0.0014353631834237914,
       0,
       0.0021530447751356872,
       0,
       -0.0014353631834237914,
       0,
       0.00035884079585594785};
  coeffs.a = coeffs_a;
  coeffs.b = coeffs_b;

  dataframe_t result;
  result.size = df->size;
  result.x = signal_filter(coeffs.size, coeffs.b, coeffs.size, coeffs.a, df->size, df->x);
  result.y = signal_filter(coeffs.size, coeffs.b, coeffs.size, coeffs.a, df->size, df->y);
  result.z = signal_filter(coeffs.size, coeffs.b, coeffs.size, coeffs.a, df->size, df->z);

  return result;
}
