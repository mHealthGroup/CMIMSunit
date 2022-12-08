#include "mims_unit.h"
#include <time.h>

dataframe_t mims_unit(dataframe_t *dataframe,
                      int dyanmic_range_low, int dyanmic_range_high,
                      int break_size, time_unit_t time_unit,
                      double noise_level, double k, double spar,
                      double cutoff_low, double cutoff_high,
                      bool allow_truncation)
{
  dataframe_t resampled_data = extrapolate(dataframe, dyanmic_range_low, dyanmic_range_high,
                                           noise_level, k, spar);

  // removed the extrapolated_data
  int sampling_rate = get_sampling_rate(&resampled_data);

  // store abnormal values separately
  int n_normal_data = 0;
  for (int i = 0; i < resampled_data.size; i++)
    if (!(resampled_data.x[i] < -150 || resampled_data.y[i] < -150 || resampled_data.z[i] < -150))
      n_normal_data += 1;

  dataframe_t normal_dataframe = {
      .size = n_normal_data,
      .x = malloc(n_normal_data * sizeof(double)),
      .y = malloc(n_normal_data * sizeof(double)),
      .z = malloc(n_normal_data * sizeof(double))};

  int *normal_rows = malloc(n_normal_data * sizeof(int));
  int normal_rows_i = 0;
  for (int i = 0; i < resampled_data.size; i++)
    if (!(resampled_data.x[i] < -150 || resampled_data.y[i] < -150 || resampled_data.z[i] < -150))
    {
      normal_rows[normal_rows_i] = i;
      normal_dataframe.x[normal_rows_i] = resampled_data.x[i];
      normal_dataframe.y[normal_rows_i] = resampled_data.y[i];
      normal_dataframe.z[normal_rows_i] = resampled_data.z[i];

      normal_rows_i += 1;
    }

  double cutoff_freq[2] = {0.2, 5.0};

  dataframe_t filtered_data = iir(&normal_dataframe, sampling_rate, cutoff_freq, 4);

  // write filtered_data back into resampled_data rows
  for (int i = 0; i < n_normal_data; i++)
  {
    resampled_data.x[normal_rows[i]] = filtered_data.x[i];
    resampled_data.y[normal_rows[i]] = filtered_data.y[i];
    resampled_data.z[normal_rows[i]] = filtered_data.z[i];
  }

  // Compute the AUC
  dataframe_t integrated_data = aggregate(&resampled_data, break_size, time_unit,
                                          true, resampled_data.timestamps[0]);

  // Truncate
  if (allow_truncation)
  {
    for (int i = 0; i < integrated_data.size; i++)
    {
      if (integrated_data.x[i] > 0 &&
          integrated_data.x[i] <= (1e-04 * parse_epoch_string(break_size, time_unit, sampling_rate)))
        integrated_data.x[i] = 0;

      if (integrated_data.y[i] > 0 &&
          integrated_data.y[i] <= (1e-04 * parse_epoch_string(break_size, time_unit, sampling_rate)))
        integrated_data.y[i] = 0;

      if (integrated_data.z[i] > 0 &&
          integrated_data.z[i] <= (1e-04 * parse_epoch_string(break_size, time_unit, sampling_rate)))
        integrated_data.z[i] = 0;
    }
  }

  integrated_data.mims_data = sum_up(&integrated_data);

  for (int i = 0; i < integrated_data.size; i++)
    if (integrated_data.x[i] < 0 || integrated_data.y[i] < 0 || integrated_data.z[i] < 0)
    {
      integrated_data.x[i] = -0.01;
      integrated_data.y[i] = -0.01;
      integrated_data.z[i] = -0.01;
      integrated_data.mims_data[i] = -0.01;
    }

  return integrated_data;
}

int main(int argc, char **argv)
{
  // int n = 108000;
  int n = 2016000;
  // int m = 720;

  int got1, got2, got3, got4;

  double *timestamps = malloc(n * sizeof(double));
  double *x = malloc(n * sizeof(double));
  double *y = malloc(n * sizeof(double));
  double *z = malloc(n * sizeof(double));
  // int *segments = malloc(m * sizeof(int));
  // FILE *t_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/mims_unit/test_2/timestamps.csv", "r");
  // FILE *x_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/mims_unit/test_2/x.csv", "r");
  // FILE *y_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/mims_unit/test_2/y.csv", "r");
  // FILE *z_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/mims_unit/test_2/z.csv", "r");

  FILE *t_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/mims_unit/aditya_sleep/timestamps.csv", "r");
  FILE *x_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/mims_unit/aditya_sleep/x.csv", "r");
  FILE *y_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/mims_unit/aditya_sleep/y.csv", "r");
  FILE *z_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/mims_unit/aditya_sleep/z.csv", "r");

  // FILE *t_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/extrapolate_input/t.csv", "r");
  // FILE *x_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/extrapolate_input/x.csv", "r");
  // FILE *y_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/extrapolate_input/y.csv", "r");
  // FILE *z_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/extrapolate_input/z.csv", "r");
  for (int i = 0; i < n; i++)
  {
    got1 = fscanf(t_file, "%lf", &timestamps[i]);
    timestamps[i] = timestamps[i] / pow(10, 9);
    got2 = fscanf(x_file, "%lf", &x[i]);
    got3 = fscanf(y_file, "%lf", &y[i]);
    got4 = fscanf(z_file, "%lf", &z[i]);
    if ((got1 + got2 + got3 + got4) != 4)
      break; // wrong number of tokens - maybe end of file
    if (x[i] != 0)
      continue;
  }
  fclose(t_file);
  fclose(x_file);
  fclose(y_file);
  fclose(z_file);

  dataframe_t input_df;
  input_df.size = n;
  input_df.timestamps = timestamps;
  input_df.x = x;
  input_df.y = y;
  input_df.z = z;

  clock_t begin = clock();
  dataframe_t mims_data = mims_unit(&input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, true);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  return 0;
}
