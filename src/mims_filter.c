#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "filter.h"
#include "helper.h"

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

// int main(int argc, char **argv)
// {
//   int sampling_rate = 100;
//   double cutoff_freq[2] = {0.2, 5.0};
//   int order = 4;
//   int n = 359997;

//   int got1, got2, got3, got4;

//   dataframe_t df;
//   df.size = n;
//   df.timestamps = malloc(n * sizeof(double));
//   df.x = malloc(n * sizeof(double));
//   df.y = malloc(n * sizeof(double));
//   df.z = malloc(n * sizeof(double));
//   FILE *iir_t_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/filter/iir_t.csv", "r");
//   FILE *iir_x_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/filter/iir_x.csv", "r");
//   FILE *iir_y_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/filter/iir_y.csv", "r");
//   FILE *iir_z_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/filter/iir_z.csv", "r");
//   for (int i = 0; i < n; i++)
//   {
//     got1 = fscanf(iir_t_file, "%lf", &df.timestamps[i]);
//     got2 = fscanf(iir_x_file, "%lf", &df.x[i]);
//     got3 = fscanf(iir_y_file, "%lf", &df.y[i]);
//     got4 = fscanf(iir_z_file, "%lf", &df.z[i]);
//     if ((got1 + got2 + got3 + got4) != 4)
//       break; // wrong number of tokens - maybe end of file
//     if (df.x[i] != 0)
//       continue;
//   }
//   fclose(iir_t_file);
//   fclose(iir_x_file);
//   fclose(iir_y_file);
//   fclose(iir_z_file);

//   dataframe_t output = iir(df, sampling_rate, cutoff_freq, order);
//   return 0;
// }
