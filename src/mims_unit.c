#include "mims_unit.h"

dataframe_t mims_unit_from_filename(char *input_filename,
                                    int8_t dyanmic_range_low, int8_t dyanmic_range_high,
                                    uint16_t break_size, time_unit_t time_unit,
                                    float noise_level, float k, float spar,
                                    float cutoff_low, float cutoff_high,
                                    uint8_t allow_truncation)
{
  dataframe_t dataframe = read_csv(input_filename);
  return custom_mims_unit(&dataframe,
                          dyanmic_range_low, dyanmic_range_high,
                          break_size, time_unit,
                          noise_level, k, spar,
                          cutoff_low, cutoff_high,
                          allow_truncation,
                          NULL, NULL);
}

dataframe_t mims_unit(dataframe_t *dataframe,
                      int8_t dyanmic_range_low, int8_t dyanmic_range_high,
                      uint16_t break_size, time_unit_t time_unit,
                      float noise_level, float k, float spar,
                      float cutoff_low, float cutoff_high,
                      uint8_t allow_truncation)
{
  return custom_mims_unit(dataframe,
                          dyanmic_range_low, dyanmic_range_high,
                          break_size, time_unit,
                          noise_level, k, spar,
                          cutoff_low, cutoff_high,
                          allow_truncation,
                          NULL, NULL);
}

dataframe_t custom_mims_unit(dataframe_t *dataframe,
                             int8_t dyanmic_range_low, int8_t dyanmic_range_high,
                             uint16_t break_size, time_unit_t time_unit,
                             float noise_level, float k, float spar,
                             float cutoff_low, float cutoff_high,
                             uint8_t allow_truncation, dataframe_t *before_df, dataframe_t *after_df)
{
  dataframe_t concatted_df;
  if (before_df)
  {
    concatted_df = concat_dataframes(before_df, dataframe);
    dataframe = &concatted_df;
  }

  if (after_df)
  {
    concatted_df = concat_dataframes(dataframe, after_df);
    dataframe = &concatted_df;
  }

  dataframe_t resampled_data = extrapolate(dataframe, dyanmic_range_low, dyanmic_range_high,
                                           noise_level, k, spar);

  // removed the extrapolated_data
  uint16_t sampling_rate = get_sampling_rate(&resampled_data);

  // store abnormal values separately
  uint32_t n_normal_data = 0;
  for (uint32_t i = 0; i < resampled_data.size; i++)
    if (!(resampled_data.x[i] < -150 || resampled_data.y[i] < -150 || resampled_data.z[i] < -150))
      n_normal_data += 1;

  dataframe_t normal_dataframe = {
      .size = n_normal_data,
      .x = malloc(n_normal_data * sizeof(double)),
      .y = malloc(n_normal_data * sizeof(double)),
      .z = malloc(n_normal_data * sizeof(double))};

  uint32_t *normal_rows = malloc(n_normal_data * sizeof(uint32_t));
  uint32_t normal_rows_i = 0;
  for (uint32_t i = 0; i < resampled_data.size; i++)
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
  for (uint32_t i = 0; i < n_normal_data; i++)
  {
    resampled_data.x[normal_rows[i]] = filtered_data.x[i];
    resampled_data.y[normal_rows[i]] = filtered_data.y[i];
    resampled_data.z[normal_rows[i]] = filtered_data.z[i];
  }

  // Compute the AUC
  dataframe_t integrated_data = aggregate(&resampled_data, break_size, time_unit,
                                          1, resampled_data.timestamps[0]);

  // Truncate
  if (allow_truncation)
  {
    for (uint32_t i = 0; i < integrated_data.size; i++)
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

  for (uint32_t i = 0; i < integrated_data.size; i++)
    if (integrated_data.x[i] < 0 || integrated_data.y[i] < 0 || integrated_data.z[i] < 0)
    {
      integrated_data.x[i] = -0.01;
      integrated_data.y[i] = -0.01;
      integrated_data.z[i] = -0.01;
      integrated_data.mims_data[i] = -0.01;
    }

  return integrated_data;
}

static dataframe_t consistency_test(char *input_filename)
{
  dataframe_t input_df = read_csv(input_filename);
  dataframe_t previous_output, current_output;
  previous_output = mims_unit(&input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1);
  for (int i = 0; i < 20; i++)
  {
    current_output = (i % 2)
                         ? mims_unit(&input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1)
                         : mims_unit_from_filename(input_filename, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1);

    for (int j = 0; j < previous_output.size; j++)
    {
      if (previous_output.mims_data[j] != current_output.mims_data[j])
      {
        printf("Failed consistency_test");
        return current_output;
      }
    }
    previous_output = current_output;
  }

  printf("Passed consistency_test");
  return current_output;
}

static void precision_test(dataframe_t *output_df, char *expected_output_filename)
{
  // Read expected outputs
  int got_r;
  int n = count_lines(expected_output_filename);
  double *r_output = malloc(n * sizeof(double));
  FILE *r_file = fopen(expected_output_filename, "r");
  for (int i = 0; i < n; i++)
  {
    got_r = fscanf(r_file, "%lf", &r_output[i]); // skip first line because it's the column name
    if (got_r != 1)
      break; // wrong number of tokens - maybe end of file
  }
  fclose(r_file);

  for (int i = 0; i < output_df->size; i++)
  {
    if (fabs(output_df->mims_data[i] - r_output[i]) > 0.0001)
    {
      printf("Failed precision test");
      return;
    }
  }

  printf("Passed precision test");
  return;
}

static void before_after_df_test(char *input_filename, char *expected_output_filename)
{
  dataframe_t input_df = read_csv(input_filename);

  dataframe_t before_df = {
      .size = 1,
      .timestamps = malloc(sizeof(double)),
      .x = malloc(sizeof(double)),
      .y = malloc(sizeof(double)),
      .z = malloc(sizeof(double))};
  dataframe_t after_df = {
      .size = 1,
      .timestamps = malloc(sizeof(double)),
      .x = malloc(sizeof(double)),
      .y = malloc(sizeof(double)),
      .z = malloc(sizeof(double))};
  before_df.timestamps[0] = input_df.timestamps[0];
  before_df.x[0] = input_df.x[0];
  before_df.y[0] = input_df.y[0];
  before_df.z[0] = input_df.z[0];

  after_df.timestamps[0] = input_df.timestamps[input_df.size - 1];
  after_df.x[0] = input_df.x[input_df.size - 1];
  after_df.y[0] = input_df.y[input_df.size - 1];
  after_df.z[0] = input_df.z[input_df.size - 1];

  dataframe_t middle_df = {
      .size = input_df.size - 2,
      .timestamps = malloc((input_df.size - 2) * sizeof(double)),
      .x = malloc((input_df.size - 2) * sizeof(double)),
      .y = malloc((input_df.size - 2) * sizeof(double)),
      .z = malloc((input_df.size - 2) * sizeof(double))};

  for (int i = 0; i < middle_df.size; i++)
  {
    middle_df.timestamps[i] = input_df.timestamps[i + 1];
    middle_df.x[i] = input_df.x[i + 1];
    middle_df.y[i] = input_df.y[i + 1];
    middle_df.z[i] = input_df.z[i + 1];
  }

  dataframe_t mims_data = custom_mims_unit(&middle_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1, &before_df, &after_df);
  precision_test(&mims_data, expected_output_filename);
}

void measure_runtime(char *input_filename)
{
#include <time.h>

  dataframe_t input_df = read_csv(input_filename);
  clock_t begin = clock();

  mims_unit(&input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  return;
}

int main(int argc, char **argv)
{
  char input_filename[66] = "/Users/arytonhoi/Kode/mhealth/cmims/data/mims_unit/test_2/raw.csv";
  char expected_output_filename[71] = "/Users/arytonhoi/Kode/mhealth/cmims/data/mims_unit/test_2/r_output.csv";

  // char filename[79] = "/Users/arytonhoi/Kode/mhealth/cmims/data/mims_unit/aditya_sleep/timestamps.csv";

  dataframe_t mims_data = consistency_test(input_filename);
  precision_test(&mims_data, expected_output_filename);
  before_after_df_test(input_filename, expected_output_filename);
  // measure_runtime(input_filename);

  return 0;
}
