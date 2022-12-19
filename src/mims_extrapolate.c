#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mims_extrapolate.h"
#include "mims_helper.h"

static values_dataframe_t extrapolate_interpolate(uint32_t n, double *oversampled_float_timestamps, double *values,
                                                  double *marker, uint16_t points_ex_n, values_dataframe_t points_ex,
                                                  uint16_t sampling_rate, float confident)
{
  uint32_t length_t_mark = 0;
  for (uint32_t i = 0; i < n; i++)
    if (fabs(marker[i]) < confident)
      length_t_mark += 1;

  uint32_t j = 0;
  uint32_t *mark_it = malloc(length_t_mark * sizeof(uint32_t));
  for (uint32_t i = 0; i < n; i++)
    if (fabs(marker[i]) < confident)
    {
      mark_it[j] = i;
      j += 1;
    }

  values_dataframe_t dat;
  if ((length_t_mark / (double)n) < 0.3)
  {
    dat.size = n;
    dat.timestamps = oversampled_float_timestamps;
    dat.values = values;
  }
  else
  {
    dat.size = length_t_mark + points_ex_n;
    dat.timestamps = malloc(dat.size * sizeof(double));
    dat.values = malloc(dat.size * sizeof(double));

    // assume dat.timestamps and points_ex.timestamps both sorted
    // mergesort into dat.timestamps
    double timestamp, value;
    uint32_t k = 0;
    j = 0;
    for (uint32_t i = 0; i < dat.size; i++)
    {
      if (
          (j < length_t_mark) &&
          ((k >= points_ex_n) || (oversampled_float_timestamps[mark_it[j]] < points_ex.timestamps[k])))
      {
        timestamp = oversampled_float_timestamps[mark_it[j]];
        value = values[mark_it[j]];
        j += 1;
      }
      else
      {
        timestamp = points_ex.timestamps[k];
        value = points_ex.values[k];
        k += 1;
      }

      dat.timestamps[i] = timestamp;
      dat.values[i] = value;
    }
  }

  free(mark_it);

  uint32_t t_interp_n = sequence_length(dat.timestamps[0], dat.timestamps[dat.size - 1], 1.0 / sampling_rate);
  double *t_interp = sequence(dat.timestamps[0], dat.timestamps[dat.size - 1], 1.0 / sampling_rate);

  // fmm spline
  double *dat_output = Spline(dat.size, dat.timestamps, dat.size, dat.values, t_interp_n, t_interp, 3);

  // Remove nan values from dat_output
  uint32_t non_nan_dat_output_n = 0;
  for (uint32_t i = 0; i < t_interp_n; i++)
  {
    if (!isnan(dat_output[i]))
      non_nan_dat_output_n += 1;
  }

  if (non_nan_dat_output_n != t_interp_n)
  {
    double *temp_dat_output = malloc(non_nan_dat_output_n * sizeof(double));
    double *temp_t_interp = malloc(non_nan_dat_output_n * sizeof(double));
    uint32_t j = 0;
    for (uint32_t i = 0; i < t_interp_n; i++)
    {
      if (!isnan(dat_output[i]))
      {
        temp_dat_output[j] = dat_output[i];
        temp_t_interp[j] = t_interp[i];
        j += 1;
      }
    }
    free(dat_output);
    free(t_interp);
    dat_output = temp_dat_output;
    t_interp_n = non_nan_dat_output_n;
  }

  // Return
  values_dataframe_t dat_df;
  dat_df.size = t_interp_n;
  free(dat.timestamps);
  free(dat.values);
  dat_df.timestamps = t_interp;
  dat_df.values = dat_output;

  return dat_df;
}

static edges_t extrapolate_edges(uint32_t n, double *marker, float confident, double sampling_rate)
{
  double *marker_diff_left = malloc(n * sizeof(double));
  double *marker_diff_right = malloc(n * sizeof(double));

  uint32_t positive_left_end_n, positive_right_start_n, negative_left_end_n, negative_right_start_n;
  positive_left_end_n = positive_right_start_n = negative_left_end_n = negative_right_start_n = 0;
  for (uint32_t i = 0; i < n; i++)
  {
    if (i < n - 1)
      marker_diff_left[i + 1] = marker_diff_right[i] = marker[i + 1] - marker[i];

    if ((marker_diff_left[i] > confident) && (marker[i] > 0))
      positive_left_end_n += 1;

    if ((marker_diff_right[i] < -confident) && (marker[i] > 0))
      positive_right_start_n += 1;

    if ((marker_diff_left[i] < -confident) && (marker[i] < 0))
      negative_left_end_n += 1;

    if ((marker_diff_right[i] > confident) && (marker[i] < 0))
      negative_right_start_n += 1;
  }

  // hills
  double *positive_left_end = malloc(positive_left_end_n * sizeof(double));
  uint32_t j;
  for (uint32_t i = 0; i < n; i++)
  {
    if ((marker_diff_left[i] > confident) && (marker[i] > 0))
    {
      positive_left_end[j] = i;
      j += 1;
    }
  }

  double *positive_right_start = malloc(positive_right_start_n * sizeof(double));
  j = 0;
  for (uint32_t i = 0; i < n; i++)
  {
    if ((marker_diff_right[i] < -confident) && (marker[i] > 0))
    {
      positive_right_start[j] = i;
      j += 1;
    }
  }

  float threshold_maxedout = sampling_rate * 5;
  if (positive_left_end_n - positive_right_start_n == 1)
  {

    if (n - positive_left_end[positive_left_end_n - 1] > threshold_maxedout) // end case > 2 second maxed out edge region
    {
      // append -1
      positive_right_start_n += 1;
      double *temp_positive_right_start = malloc(positive_right_start_n * sizeof(double));
      memcpy(temp_positive_right_start, positive_right_start, (positive_right_start_n - 1) * sizeof(double));
      free(positive_right_start);
      temp_positive_right_start[positive_right_start_n - 1] = -1;
      positive_right_start = temp_positive_right_start;
    }
    else
    {
      // < 2 second maxed out edge region, do nothing
      if (positive_left_end_n == 1)
        positive_left_end_n = 0;
      else // "remove" last element
        positive_left_end_n -= 1;
    }
  }
  else if (positive_left_end_n - positive_right_start_n == -1)
  { // start case > 2 second maxed out edge region
    if (positive_right_start[0] > threshold_maxedout)
    {
      // prepend -1 to positive_left_end
      positive_left_end_n += 1;
      double *temp_positive_left_end = malloc(positive_left_end_n * sizeof(double));
      memcpy(
          temp_positive_left_end + 1,
          positive_left_end,
          (positive_left_end_n - 1) * sizeof(double));
      free(positive_left_end);
      temp_positive_left_end[0] = -1;
      positive_left_end = temp_positive_left_end;
    }
    else
    {
      // < 2 second maxed out edge region, do nothing
      if (positive_right_start_n == 1)
        positive_right_start_n = 0;
      else
      {
        // "remove" first element
        positive_right_start_n -= 1;
        positive_right_start += 1;
      }
    }
  }

  // Why would array of indices have nan? These two lines seem useless
  // positive_left_end[~np.isnan(positive_left_end)]
  // positive_right_start[~np.isnan(positive_right_start)]

  uint8_t positive_right_start_lessthan_right_end = 0;
  for (uint32_t i = 0; i < min(positive_left_end_n, positive_right_start_n); i++)
    if (positive_right_start[i] < positive_left_end[i])
    {
      positive_right_start_lessthan_right_end = 1;
      break;
    }

  if (positive_right_start_n > 1 && positive_right_start_lessthan_right_end)
  {
    positive_left_end_n -= 1;
    positive_right_start_n -= 1; // "Remove first element"
    positive_right_start += 1;
  }

  // Shouldn't need this. Was only used in python because pandas.df wanted
  // vertical arrays or something
  // positive_edges = np.column_stack((positive_left_end, positive_right_start))

  // valleys
  double *negative_left_end = malloc(negative_left_end_n * sizeof(double));
  j = 0;
  for (uint32_t i = 0; i < n; i++)
  {
    if ((marker_diff_left[i] < -confident) && (marker[i] < 0))
    {
      negative_left_end[j] = i;
      j += 1;
    }
  }

  double *negative_right_start = malloc(negative_right_start_n * sizeof(double));
  j = 0;
  for (uint32_t i = 0; i < n; i++)
  {
    if ((marker_diff_right[i] > confident) && (marker[i] < 0))
    {
      negative_right_start[j] = i;
      j += 1;
    }
  }

  if (negative_left_end_n - negative_right_start_n == 1)
  {
    if (n - negative_left_end[negative_left_end_n - 1] > threshold_maxedout) // end case > 2 second maxed out edge region
    {
      // append -1
      negative_right_start_n += 1;
      double *temp_negative_right_start = malloc(negative_right_start_n * sizeof(double));
      memcpy(temp_negative_right_start, negative_right_start, (negative_right_start_n - 1) * sizeof(double));
      free(negative_right_start);
      temp_negative_right_start[negative_right_start_n - 1] = -1;
      negative_right_start = temp_negative_right_start;
    }
    else
    {
      // < 2 second maxed out edge region, do nothing
      if (negative_left_end_n == 1)
        negative_left_end_n = 0;
      else // "remove" last element
        negative_left_end_n -= 1;
    }
  }
  else if (negative_left_end_n - negative_right_start_n == -1)
  { // start case > 5 second maxed out edge region
    if (negative_right_start[0] > threshold_maxedout)
    {
      // prepend -1 to negative_left_end
      negative_left_end_n += 1;
      double *temp_negative_left_end = malloc(negative_left_end_n * sizeof(double));
      memcpy(
          temp_negative_left_end + 1,
          negative_left_end,
          (negative_left_end_n - 1) * sizeof(double));
      free(negative_left_end);
      temp_negative_left_end[0] = -1;
      negative_left_end = temp_negative_left_end;
    }
    else
    {
      // < 2 second maxed out edge region, do nothing
      if (negative_right_start_n == 1)
        negative_right_start_n = 0;
      else
      {
        // "remove" first element
        negative_right_start_n -= 1;
        negative_right_start += 1;
      }
    }
  }

  // Why would array of indices have nan? These two lines seem useless
  // negative_left_end[~np.isnan(negative_left_end)]
  // negative_right_start[~np.isnan(negative_right_start)]

  uint8_t negative_right_start_lessthan_right_end = 0;
  for (uint32_t i = 0; i < min(negative_left_end_n, negative_right_start_n); i++)
    if (negative_right_start[i] < negative_left_end[i])
    {
      negative_right_start_lessthan_right_end = 1;
      break;
    }

  if (negative_right_start_n > 1 && negative_right_start_lessthan_right_end)
  {
    negative_left_end_n -= 1;    // "remove last element"
    negative_right_start_n -= 1; // "Remove first element"
    negative_right_start += 1;
  }

  // Shouldn't need this. Was only used in python because pandas.df wanted
  // vertical arrays or something
  // negative_edges = np.column_stack((negative_left_end, negative_right_start))

  free(marker_diff_left);
  free(marker_diff_right);

  edges_t edges;
  edges.n_left = positive_left_end_n + negative_left_end_n;
  edges.left_end = malloc(edges.n_left * sizeof(uint32_t));
  for (uint16_t i = 0; i < edges.n_left; i++)
  {
    if (i < positive_left_end_n)
      edges.left_end[i] = positive_left_end[i];
    else
      edges.left_end[i] = negative_left_end[i - positive_left_end_n];
  }
  free(positive_left_end);
  free(negative_left_end);

  edges.n_right = positive_right_start_n + negative_right_start_n;
  edges.right_start = malloc(edges.n_right * sizeof(uint32_t));
  for (uint16_t i = 0; i < edges.n_right; i++)
  {
    if (i < positive_right_start_n)
      edges.right_start[i] = positive_right_start[i];
    else
      edges.right_start[i] = negative_right_start[i - positive_right_start_n];
  }

  free(positive_right_start);
  free(negative_right_start);

  return edges;
}

static edges_t extrapolate_neighbor(uint32_t n, double *marker, double sampling_rate, float k, float confident)
{
  uint16_t n_neighbor = (uint16_t)(k * sampling_rate);
  edges_t edges = extrapolate_edges(n, marker, confident, sampling_rate);

  if (edges.n_left > 0)
  {
    edges.left_start = malloc(edges.n_left * sizeof(uint32_t));
    for (uint16_t i = 0; i < edges.n_left; i++)
      edges.left_start[i] = (edges.left_end[i] == -1) ? -1 : max(edges.left_end[i] - n_neighbor + 1, 1);
  }

  if (edges.n_right > 0)
  {
    edges.right_end = malloc(edges.n_right * sizeof(uint32_t));
    for (uint16_t i = 0; i < edges.n_right; i++)
      edges.right_end[i] = (edges.right_start[i] == -1) ? -1 : min(edges.right_start[i] + n_neighbor - 1, n);
  }

  return edges;
}

static smooth_spline_model_t fit_weighted(uint32_t n, double *oversampled_float_timestamps, double *values,
                                          double *marker, uint32_t start, uint32_t end, float spar,
                                          uint16_t sampling_rate, float k)
{
  uint16_t n_over = k * sampling_rate;

  uint32_t n_sub = end - start + 1;
  double *sub_timestamps = malloc(n_sub * sizeof(double));
  double *sub_values = malloc(n_sub * sizeof(double));
  double *sub_weights = malloc(n_sub * sizeof(double));
  for (uint32_t i = 0; i < n_sub; i++)
  {
    sub_timestamps[i] = oversampled_float_timestamps[start + i];
    sub_values[i] = values[start + i];
    sub_weights[i] = 1.0 - marker[start + i];
  }

  approx_output_t sp = approx(n_sub, sub_timestamps, sub_values, n_over);
  approx_output_t weights = approx(n_sub, sub_timestamps, sub_weights, n_over);

  double *over_t = sp.x;
  double *over_value = sp.y;
  double *weight = weights.y;

  return SmSplineCoef(n_over, over_t, over_value, n_over, weight, spar);
}

static values_dataframe_t extrapolate_fit(uint32_t n, double *oversampled_float_timestamps, double *values,
                                          edges_t neighbors, double *marker, double spar, uint16_t sampling_rate,
                                          double k)
{
  // assuming neighbors.n_left == neighbors.n_right
  double *middle_ts = malloc(neighbors.n_left * sizeof(double));
  double *point_exs = malloc(neighbors.n_left * sizeof(double));

  smooth_spline_model_t fitted_left;
  smooth_spline_model_t fitted_right;
  double *left_ex, *right_ex, *middle_t, point_ex, start_time, left_end, right_start;
  middle_t = malloc(sizeof(double));
  for (uint16_t i = 0; i < neighbors.n_left; i++)
  {
    fitted_left = fit_weighted(n, oversampled_float_timestamps, values, marker,
                               neighbors.left_start[i], neighbors.left_end[i],
                               spar, sampling_rate, k);
    fitted_right = fit_weighted(n, oversampled_float_timestamps, values, marker,
                                neighbors.right_start[i], neighbors.right_end[i],
                                spar, sampling_rate, k);

    start_time = oversampled_float_timestamps[neighbors.left_end[i]];
    left_end = 0;
    right_start = oversampled_float_timestamps[neighbors.right_start[i]] - start_time;

    middle_t[0] = ((double)(left_end + right_start) / 2.0) + start_time;
    left_ex = predict_smooth_spline(fitted_left, middle_t, 1, 0);
    right_ex = predict_smooth_spline(fitted_right, middle_t, 1, 0);
    point_ex = (left_ex[0] + right_ex[0]) / 2.0;

    middle_ts[i] = middle_t[0];
    point_exs[i] = point_ex;
  }

  values_dataframe_t point_exs_df = {
      .size = neighbors.n_left,
      .timestamps = middle_ts,
      .values = point_exs};

  return point_exs_df;
}

static double optimize_gamma(double value)
{
  double start = 0.5;
  double stop = 0.001;
  double step = -0.001;

  uint32_t n = sequence_length(0.5, 0.001, -0.001);
  double *arr = sequence(0.5, 0.001, -0.001);

  double ii, current;
  double result = 0;
  double previous = 1.0;
  double previous_ii = 0.0;
  for (uint32_t i = 0; i < n; i++)
  {
    ii = arr[i];
    current = pgamma(value, ii, 1.0);
    if (previous < 0.95 && current >= 0.95)
    {
      result = (fabs(0.95 - previous) > fabs(current - 0.95)) ? ii : previous_ii;
      break;
    }
    previous = current;
    previous_ii = ii;
  }
  return result;
}

static double *mark_gamma(uint32_t n, double *timestamps, double *values, int8_t range_low, int8_t range_high,
                          float noise_sd)
{
  double *marker = calloc(n, sizeof(double));

  // model of 3sd and shape para with confident probability at 0.95
  noise_sd += 0.00001;
  double shape = optimize_gamma(3.0 * noise_sd);

  // mark using gamma distribution
  for (uint32_t i = 0; i < n; i++)
  {
    if (values[i] >= 0)
      marker[i] = pgamma(values[i] - (range_high - 5.0 * noise_sd), shape, 1.0);
    else
      marker[i] = -1.0 * pgamma(-1.0 * values[i] + (range_low + 5.0 * noise_sd), shape, 1.0);
  }

  return marker;
}

static values_dataframe_t extrapolate_single_col(uint32_t n, double *timestamps, double *values,
                                                 int8_t r_low, int8_t r_high, float noise_level, float k, float spar)
{

  uint32_t oversampled_float_timestamps_n = sequence_length(timestamps[0], timestamps[n - 1], 0.01);
  double *oversampled_float_timestamps = sequence(timestamps[0], timestamps[n - 1], 0.01);

  // method = natural
  // assume all oversampled_float_timestamps unique (python code)
  // _, indices = np.unique(oversampled_float_timestamps, return_index=True)
  // values = dat_over[indices]

  double *dat_over = Spline(n, timestamps, n, values, oversampled_float_timestamps_n, oversampled_float_timestamps, 2);
  values = dat_over;

  // mark maxed out region using gamma distribution or threshold
  double *marker = mark_gamma(oversampled_float_timestamps_n, oversampled_float_timestamps, values, r_low, r_high, noise_level);

  // mark neighbors
  edges_t neighbors = extrapolate_neighbor(oversampled_float_timestamps_n, marker, 100, k, 0.5);

  // fit local spline regression
  values_dataframe_t points_ex = extrapolate_fit(oversampled_float_timestamps_n, oversampled_float_timestamps, values,
                                                 neighbors, marker, spar, 100, k);

  // interpolate with the original
  values_dataframe_t dat_interp = extrapolate_interpolate(
      oversampled_float_timestamps_n,
      oversampled_float_timestamps,
      values,
      marker,
      points_ex.size,
      points_ex,
      100,
      0.5);

  return dat_interp;
}

dataframe_t extrapolate(dataframe_t *df, int8_t r_low, int8_t r_high, float noise_level, float k, float spar)
{
  dataframe_t result;
  values_dataframe_t x_col = extrapolate_single_col(df->size, df->timestamps, df->x, r_low, r_high, noise_level, k, spar);
  values_dataframe_t y_col = extrapolate_single_col(df->size, df->timestamps, df->y, r_low, r_high, noise_level, k, spar);
  values_dataframe_t z_col = extrapolate_single_col(df->size, df->timestamps, df->z, r_low, r_high, noise_level, k, spar);

  result.size = x_col.size;
  result.timestamps = x_col.timestamps;
  result.x = x_col.values;
  result.y = y_col.values;
  result.z = z_col.values;

  return result;
}
