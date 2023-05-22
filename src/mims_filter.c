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
dataframe_t *iir(dataframe_t *df, const uint16_t sampling_rate, double *cutoff_freq, const uint8_t order)
{
  // uint8_t cutoff_freq_n = 2;
  // float nyquist = (float)sampling_rate / 2.0;
  // for (uint8_t i = 0; i < cutoff_freq_n; i++)
  //   cutoff_freq[i] = cutoff_freq[i] / nyquist;

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

  transfer_t *coeffs = malloc(sizeof(transfer_t));
  coeffs->size = 9;
  coeffs->a = malloc(coeffs->size * sizeof(double));
  coeffs->b = malloc(coeffs->size * sizeof(double));

  coeffs->a[0] = 1.0;
  coeffs->a[1] = -7.1989557602839946;
  coeffs->a[2] = 22.71024005695034;
  coeffs->a[3] = -41.012138781845756;
  coeffs->a[4] = 46.37860384275109;
  coeffs->a[5] = -33.634144664989293;
  coeffs->a[6] = 15.276672513261936;
  coeffs->a[7] = -3.9733824787545187;
  coeffs->a[8] = 0.4531052730785417;

  coeffs->b[0] = 0.00035884079585594785;
  coeffs->b[1] = 0;
  coeffs->b[2] = -0.0014353631834237914;
  coeffs->b[3] = 0;
  coeffs->b[4] = 0.0021530447751356872;
  coeffs->b[5] = 0;
  coeffs->b[6] = -0.0014353631834237914;
  coeffs->b[7] = 0;
  coeffs->b[8] = 0.00035884079585594785;

  dataframe_t *result = create_dataframe(
      df->size,
      NULL,
      signal_filter(coeffs->size, coeffs->b, coeffs->size, coeffs->a, df->size, df->x),
      signal_filter(coeffs->size, coeffs->b, coeffs->size, coeffs->a, df->size, df->y),
      signal_filter(coeffs->size, coeffs->b, coeffs->size, coeffs->a, df->size, df->z),
      0, NULL, NULL);

  free(coeffs->a);
  free(coeffs->b);
  free(coeffs);

  return result;
}
