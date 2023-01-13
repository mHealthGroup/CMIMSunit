#include "mims_unit.h"
#include "config.h"
#include <time.h>
#include "mims_helper.h"

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
