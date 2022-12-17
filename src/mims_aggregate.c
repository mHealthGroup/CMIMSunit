#include "mims_aggregate.h"

// note: run integrate_for_mims on each column (x, y, z) separately
void integrate_for_mims(double *result, double *new_timestamps, int segments_n, int *segments,
                        int n, double *float_timestamps, double *values, int n_threshold,
                        uint8_t rectify)
{
  int segment_start_i, segment_end_i, segment_length, max_values;
  int result_i = 0;
  double auc_value;
  for (int segment_i = 0; segment_i < segments_n; segment_i++)
  {
    segment_start_i = segments[segment_i];
    segment_end_i = (segment_i == (segments_n - 1)) ? n : segments[segment_i + 1];
    segment_length = segment_end_i - segment_start_i;

    if (segment_length >= (0.9 * n_threshold))
    {
      if (rectify) // absolute values for values > constant
        for (int j = segment_start_i; j < segment_end_i; j++)
          values[j] = (values[j] > -150) ? fabs(values[j]) : -200;

      // select different methods for integration
      auc_value = catools_trapz(
          segment_length,
          float_timestamps + segment_start_i,
          values + segment_start_i);

      max_values = 16 * n_threshold;
    }
    else
    {
      auc_value = -1;
      max_values = 0;
    }

    // flag extra large (abnormal) values
    if ((auc_value >= max_values) || (auc_value < 0))
      auc_value = -1;

    // write to result array
    result[result_i] = auc_value;
    new_timestamps[result_i] = float_timestamps[segment_start_i];
    result_i += 1;
  }

  return;
}

dataframe_t aggregate(dataframe_t *dataframe, int break_size, time_unit_t time_unit,
                      uint8_t rectify, double start_time)
{
  // parse input argument epoch
  segment_data(dataframe, break_size, time_unit, start_time);

  // get the number of samples in each epoch
  int sampling_rate = get_sampling_rate(dataframe);
  int n_threshold = parse_epoch_string(break_size, time_unit, sampling_rate);

  double *x_results = malloc(dataframe->n_segments * sizeof(double));
  double *new_timestamps = malloc(dataframe->n_segments * sizeof(double));
  integrate_for_mims(x_results, new_timestamps, dataframe->n_segments, dataframe->segments,
                     dataframe->size, dataframe->timestamps, dataframe->x, n_threshold, 1);

  double *y_results = malloc(dataframe->n_segments * sizeof(double));
  integrate_for_mims(y_results, new_timestamps, dataframe->n_segments, dataframe->segments,
                     dataframe->size, dataframe->timestamps, dataframe->y, n_threshold, 1);

  double *z_results = malloc(dataframe->n_segments * sizeof(double));
  integrate_for_mims(z_results, new_timestamps, dataframe->n_segments, dataframe->segments,
                     dataframe->size, dataframe->timestamps, dataframe->z, n_threshold, 1);

  dataframe_t new_dataframe;
  new_dataframe.size = dataframe->n_segments;
  new_dataframe.timestamps = new_timestamps;
  new_dataframe.x = x_results;
  new_dataframe.y = y_results;
  new_dataframe.z = z_results;

  return new_dataframe;
}