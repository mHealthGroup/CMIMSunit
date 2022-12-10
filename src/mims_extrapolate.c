#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "mims_extrapolate.h"
#include "mims_helper.h"

static values_dataframe_t extrapolate_interpolate(int n, double *oversampled_float_timestamps, double *values,
                                                  double *marker, int points_ex_n, values_dataframe_t points_ex,
                                                  int sampling_rate, double confident)
{
  int length_t_mark = 0;
  for (int i = 0; i < n; i++)
    if (fabs(marker[i]) < confident)
      length_t_mark += 1;

  int j = 0;
  int *mark_it = malloc(length_t_mark * sizeof(int));
  for (int i = 0; i < n; i++)
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
    int k = 0;
    j = 0;
    for (int i = 0; i < dat.size; i++)
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

  int t_interp_n = sequence_length(dat.timestamps[0], dat.timestamps[dat.size - 1], 1.0 / sampling_rate);
  double *t_interp = sequence(dat.timestamps[0], dat.timestamps[dat.size - 1], 1.0 / sampling_rate);

  // DELETE ---------------------
  int test_n, got1;
  FILE *t_file, *v_file;

  // test_n = 359997;
  // t_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/extra_inter_spline_output/t.csv", "r");
  // t_interp_n = test_n;
  // free(t_interp);
  // t_interp = malloc(test_n * sizeof(double));
  // for (int i = 0; i < test_n; i++)
  // {
  //   got1 = fscanf(t_file, "%lf", &t_interp[i]);
  //   if (got1 != 1)
  //     break; // wrong number of tokens - maybe end of file
  //   if (t_interp[i] != 0)
  //     continue;
  // }
  // fclose(t_file);
  // ---------------------

  // fmm spline
  double *dat_output = Spline(dat.size, dat.timestamps, dat.size, dat.values, t_interp_n, t_interp, 3);

  // DELETE ---------------------

  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/c/extra_inter_spline_output/v.csv", "w+");
  // for (int i = 0; i < t_interp_n; i++)
  // {
  //   fprintf(v_file, "%lf\n", dat_output[i]);
  // }
  // fclose(v_file);

  test_n = 359997;
  // free(dat_output);
  // dat_output = malloc(test_n * sizeof(double));
  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/extra_inter_spline_output/v.csv", "r");
  // for (int i = 0; i < test_n; i++)
  // {
  //   got1 = fscanf(v_file, "%lf", &dat_output[i]);
  //   if (got1 != 1)
  //     break; // wrong number of tokens - maybe end of file
  //   if (dat_output[i] != 0)
  //     continue;
  // }
  // fclose(v_file);
  // ---------------------

  // Remove nan values from dat_output
  int non_nan_dat_output_n = 0;
  for (int i = 0; i < t_interp_n; i++)
  {
    if (!isnan(dat_output[i]))
      non_nan_dat_output_n += 1;
  }

  if (non_nan_dat_output_n != t_interp_n)
  {
    double *temp_dat_output = malloc(non_nan_dat_output_n * sizeof(double));
    double *temp_t_interp = malloc(non_nan_dat_output_n * sizeof(double));
    int j = 0;
    for (int i = 0; i < t_interp_n; i++)
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
  // for (int i = 0; i < t_interp_n; i++)
  //   t_interp[i] = t_interp[i] * pow(10, 9);
  free(dat.timestamps);
  free(dat.values);
  dat_df.timestamps = t_interp;
  dat_df.values = dat_output;

  return dat_df;
}

static edges_t extrapolate_edges(int n, double *marker, double confident, double sampling_rate)
{
  double *marker_diff_left = malloc(n * sizeof(double));
  double *marker_diff_right = malloc(n * sizeof(double));

  int positive_left_end_n, positive_right_start_n, negative_left_end_n, negative_right_start_n;
  positive_left_end_n = positive_right_start_n = negative_left_end_n = negative_right_start_n = 0;
  for (int i = 0; i < n; i++)
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
  int j;
  for (int i = 0; i < n; i++)
  {
    if ((marker_diff_left[i] > confident) && (marker[i] > 0))
    {
      positive_left_end[j] = i;
      j += 1;
    }
  }

  double *positive_right_start = malloc(positive_right_start_n * sizeof(double));
  j = 0;
  for (int i = 0; i < n; i++)
  {
    if ((marker_diff_right[i] < -confident) && (marker[i] > 0))
    {
      positive_right_start[j] = i;
      j += 1;
    }
  }

  double threshold_maxedout = sampling_rate * 5;
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

  bool positive_right_start_lessthan_right_end = false;
  for (int i = 0; i < min(positive_left_end_n, positive_right_start_n); i++)
    if (positive_right_start[i] < positive_left_end[i])
    {
      positive_right_start_lessthan_right_end = true;
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
  for (int i = 0; i < n; i++)
  {
    if ((marker_diff_left[i] < -confident) && (marker[i] < 0))
    {
      negative_left_end[j] = i;
      j += 1;
    }
  }

  double *negative_right_start = malloc(negative_right_start_n * sizeof(double));
  j = 0;
  for (int i = 0; i < n; i++)
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

  bool negative_right_start_lessthan_right_end = false;
  for (int i = 0; i < min(negative_left_end_n, negative_right_start_n); i++)
    if (negative_right_start[i] < negative_left_end[i])
    {
      negative_right_start_lessthan_right_end = true;
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
  edges.left_end = malloc(edges.n_left * sizeof(int));
  for (int i = 0; i < edges.n_left; i++)
  {
    if (i < positive_left_end_n)
      edges.left_end[i] = positive_left_end[i];
    else
      edges.left_end[i] = negative_left_end[i - positive_left_end_n];
  }
  free(positive_left_end);
  free(negative_left_end);

  edges.n_right = positive_right_start_n + negative_right_start_n;
  edges.right_start = malloc(edges.n_right * sizeof(int));
  for (int i = 0; i < edges.n_right; i++)
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

static edges_t extrapolate_neighbor(int n, double *marker, double sampling_rate, double k, double confident)
{
  int n_neighbor = (int)(k * sampling_rate);
  edges_t edges = extrapolate_edges(n, marker, confident, sampling_rate);

  if (edges.n_left > 0)
  {
    edges.left_start = malloc(edges.n_left * sizeof(int));
    for (int i = 0; i < edges.n_left; i++)
      edges.left_start[i] = (edges.left_end[i] == -1) ? -1 : max(edges.left_end[i] - n_neighbor + 1, 1);
  }

  if (edges.n_right > 0)
  {
    edges.right_end = malloc(edges.n_right * sizeof(int));
    for (int i = 0; i < edges.n_right; i++)
      edges.right_end[i] = (edges.right_start[i] == -1) ? -1 : min(edges.right_start[i] + n_neighbor - 1, n);
  }

  return edges;
}

static smooth_spline_model_t fit_weighted(int n, double *oversampled_float_timestamps, double *values,
                                          double *marker, int start, int end, double spar,
                                          int sampling_rate, double k)
{
  int n_over = k * sampling_rate;

  int n_sub = end - start + 1;
  double *sub_timestamps = malloc(n_sub * sizeof(double));
  double *sub_values = malloc(n_sub * sizeof(double));
  double *sub_weights = malloc(n_sub * sizeof(double));
  for (int i = 0; i < n_sub; i++)
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

static values_dataframe_t extrapolate_fit(int n, double *oversampled_float_timestamps, double *values, edges_t neighbors,
                                          double *marker, double spar, int sampling_rate, double k)
{
  // assuming neighbors.n_left == neighbors.n_right
  double *middle_ts = malloc(neighbors.n_left * sizeof(double));
  double *point_exs = malloc(neighbors.n_left * sizeof(double));

  smooth_spline_model_t fitted_left;
  smooth_spline_model_t fitted_right;
  double *left_ex, *right_ex, *middle_t, point_ex, start_time, left_end, right_start;
  middle_t = malloc(sizeof(double));
  for (int i = 0; i < neighbors.n_left; i++)
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
  // int n = (int)((0.001 - 0.5) / -0.001) + 1;
  int n = sequence_length(0.5, 0.001, -0.001);
  double *arr = sequence(0.5, 0.001, -0.001);

  double ii, current;
  double result = 0;
  double previous = 1.0;
  double previous_ii = 0.0;
  for (int i = 0; i < n; i++)
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

static double *mark_gamma(int n, double *timestamps, double *values, double range_low, double range_high,
                          double noise_sd)
{
  double *marker = calloc(n, sizeof(double));

  // model of 3sd and shape para with confident probability at 0.95
  noise_sd += 0.00001;
  double shape = optimize_gamma(3.0 * noise_sd);

  // mark using gamma distribution
  for (int i = 0; i < n; i++)
  {
    if (values[i] >= 0)
      marker[i] = pgamma(values[i] - (range_high - 5.0 * noise_sd), shape, 1.0);
    else
      marker[i] = -1.0 * pgamma(-1.0 * values[i] + (range_low + 5.0 * noise_sd), shape, 1.0);
  }

  return marker;
}

static values_dataframe_t extrapolate_single_col(int n, double *timestamps, double *values,
                                                 double r_low, double r_high, double noise_level, double k, double spar)
{

  int oversampled_float_timestamps_n = sequence_length(timestamps[0], timestamps[n - 1], 0.01);
  double *oversampled_float_timestamps = sequence(timestamps[0], timestamps[n - 1], 0.01);

  // method = natural
  // assume all oversampled_float_timestamps unique (python code)
  // _, indices = np.unique(oversampled_float_timestamps, return_index=True)
  // values = dat_over[indices]

  // DELETE ---------------------
  int test_n, got1, got2;
  FILE *t_file, *v_file;
  // ---------------------

  // DELETE ---------------------
  // t_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/c/spline_input/t.csv", "w+");
  // for (int i = 0; i < oversampled_float_timestamps_n; i++)
  // {
  //   fprintf(t_file, "%lf\n", oversampled_float_timestamps[i]);
  // }
  // fclose(t_file);
  test_n = 359997;

  // t_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/spline_input/over_t.csv", "r");
  // oversampled_float_timestamps_n = test_n;
  // free(oversampled_float_timestamps);
  // oversampled_float_timestamps = malloc(test_n * sizeof(double));
  // for (int i = 0; i < test_n; i++)
  // {
  //   got1 = fscanf(t_file, "%lf", &oversampled_float_timestamps[i]);
  //   if (got1 != 1)
  //     break; // wrong number of tokens - maybe end of file
  //   if (oversampled_float_timestamps[i] != 0)
  //     continue;
  // }
  // fclose(t_file);
  // ---------------------

  double *dat_over = Spline(n, timestamps, n, values, oversampled_float_timestamps_n, oversampled_float_timestamps, 2);
  values = dat_over;

  // DELETE ---------------------

  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/c/spline_output/v.csv", "w+");
  // for (int i = 0; i < oversampled_float_timestamps_n; i++)
  // {
  //   fprintf(v_file, "%lf\n", values[i]);
  // }
  // fclose(v_file);

  test_n = 359997;
  // free(values);
  // values = malloc(test_n * sizeof(double));
  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/spline_output/v.csv", "r");
  // for (int i = 0; i < test_n; i++)
  // {
  //   got1 = fscanf(v_file, "%lf", &values[i]);
  //   if (got1 != 1)
  //     break; // wrong number of tokens - maybe end of file
  //   if (values[i] != 0)
  //     continue;
  // }
  // fclose(v_file);
  // ---------------------

  // mark maxed out region using gamma distribution or threshold
  double *marker = mark_gamma(oversampled_float_timestamps_n, oversampled_float_timestamps, values, r_low, r_high, noise_level);
  // DELETE ---------------------
  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/c/mark_output/v.csv", "w+");
  // for (int i = 0; i < oversampled_float_timestamps_n; i++)
  // {
  //   fprintf(v_file, "%lf\n", marker[i]);
  // }
  // fclose(v_file);

  test_n = 359997;
  // free(marker);
  // marker = malloc(test_n * sizeof(double));
  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/mark_output/v.csv", "r");
  // for (int i = 0; i < test_n; i++)
  // {
  //   got1 = fscanf(v_file, "%lf", &marker[i]);
  //   if (got1 != 1)
  //     break; // wrong number of tokens - maybe end of file
  //   if (marker[i] != 0)
  //     continue;
  // }
  // fclose(v_file);
  // ---------------------

  // mark neighbors
  edges_t neighbors = extrapolate_neighbor(oversampled_float_timestamps_n, marker, 100, k, 0.5);

  // fit local spline regression
  values_dataframe_t points_ex = extrapolate_fit(oversampled_float_timestamps_n, oversampled_float_timestamps, values,
                                                 neighbors, marker, spar, 100, k);
  // DELETE ---------------------
  test_n = 0;

  // t_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/c/points_ex_output/t.csv", "w+");
  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/c/points_ex_output/v.csv", "w+");
  // for (int i = 0; i < points_ex.size; i++)
  // {
  //   fprintf(t_file, "%lf\n", points_ex.timestamps[i]);
  //   fprintf(v_file, "%lf\n", points_ex.values[i]);
  // }
  // fclose(t_file);
  // fclose(v_file);

  // free(points_ex.timestamps);
  // free(points_ex.values);
  // points_ex.timestamps = malloc(test_n * sizeof(double));
  // points_ex.values = malloc(test_n * sizeof(double));
  // t_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/points_ex_output/t.csv", "r");
  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/points_ex_output/v.csv", "r");
  // for (int i = 0; i < test_n; i++)
  // {
  //   got1 = fscanf(t_file, "%lf", &points_ex.timestamps[i]);
  //   got2 = fscanf(v_file, "%lf", &points_ex.values[i]);
  //   if ((got1 + got2) != 2)
  //     break; // wrong number of tokens - maybe end of file
  //   if (points_ex.values[i] != 0)
  //     continue;
  // }
  // fclose(t_file);
  // fclose(v_file);
  // ---------------------

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
  // DELETE ---------------------

  // t_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/c/extra_inter_output/t.csv", "w+");
  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/c/extra_inter_output/v.csv", "w+");
  // for (int i = 0; i < dat_interp.size; i++)
  // {
  //   fprintf(t_file, "%lf\n", dat_interp.timestamps[i]);
  //   fprintf(v_file, "%lf\n", dat_interp.values[i]);
  // }
  // fclose(t_file);
  // fclose(v_file);

  test_n = 359997;
  // free(dat_interp.timestamps);
  // free(dat_interp.values);
  // dat_interp.timestamps = malloc(test_n * sizeof(double));
  // dat_interp.values = malloc(test_n * sizeof(double));
  // t_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/extra_inter_output/t.csv", "r");
  // v_file = fopen("/Users/arytonhoi/Kode/mhealth/pymims/test-data/cmims/r/extra_inter_output/v.csv", "r");
  // for (int i = 0; i < test_n; i++)
  // {
  //   got1 = fscanf(t_file, "%lf", &dat_interp.timestamps[i]);
  //   got2 = fscanf(v_file, "%lf", &dat_interp.values[i]);
  //   if ((got1 + got2) != 2)
  //     break; // wrong number of tokens - maybe end of file
  //   if (dat_interp.values[i] != 0)
  //     continue;
  // }
  // fclose(t_file);
  // fclose(v_file);
  // ---------------------

  return dat_interp;
}

dataframe_t extrapolate(dataframe_t *df, double r_low, double r_high, double noise_level, double k, double spar)
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

void run_subtests()
{
  int num_tests = 3;
  int m;
  double *test_array, *temp_test_array;
  for (int test_i = 0; test_i < num_tests; test_i++)
  {
    int m = 10;
    test_array = malloc(m * sizeof(double));
    for (int i = 0; i < 10; i++)
      test_array[i] = i;

    switch (test_i)
    {
    case 0: // Test appending -1
      m += 1;
      temp_test_array = malloc(m * sizeof(double));
      memcpy(temp_test_array, test_array, (m - 1) * sizeof(double));
      temp_test_array[m - 1] = -1;
      free(test_array);
      test_array = temp_test_array;
      free(test_array);
      break;

    case 1: // Test preprending -1
      m += 1;
      temp_test_array = malloc(m * sizeof(double));
      memcpy(
          temp_test_array + 1,
          test_array,
          (m - 1) * sizeof(double));
      temp_test_array[0] = -1;
      free(test_array);
      test_array = temp_test_array;
      free(test_array);
      break;

    case 2: // Test removing first element
      m -= 1;
      test_array += 1;
      break;

    default:
      break;
    }
  }
}

// int main(int argc, char **argv)
// {
//   int n = 359997;
//   int n2 = 108000;
//   // double confident = 0.5;
//   // int sampling_rate = 100;
//   // double k = 0.05;
//   // double spar = 0.6;

//   // double *edges_marker = malloc(n * sizeof(double));

//   // double *interpolate_marker = malloc(n * sizeof(double));
//   // double *interpolate_oversampled_float_timestamps = malloc(n * sizeof(double));
//   // double *interpolate_values = malloc(n * sizeof(double));

//   // double *fit_weighted_over_t = malloc(n * sizeof(double));
//   // double *fit_weighted_values = malloc(n * sizeof(double));
//   // double *fit_weighted_marker = malloc(n * sizeof(double));

//   // double *mark_gamma_t = malloc(n * sizeof(double));
//   // double *mark_gamma_values = malloc(n * sizeof(double));

//   double *single_col_t = malloc(n2 * sizeof(double));
//   double *single_col_values = malloc(n2 * sizeof(double));

//   int got1, got2, got3;

//   // FILE *edges_marker_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/edges_marker.csv", "r");
//   // FILE *interpolate_marker_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/interpolate_marker.csv", "r");
//   // for (int i = 0; i < n; i++)
//   // {
//   //   got1 = fscanf(edges_marker_file, "%lf", &edges_marker[i]);
//   //   got2 = fscanf(interpolate_marker_file, "%lf", &interpolate_marker[i]);
//   //   if ((got1 + got2) != 2)
//   //     break; // wrong number of tokens - maybe end of file
//   //   if (edges_marker[i] != 0)
//   //     continue;
//   // }
//   // fclose(edges_marker_file);
//   // fclose(interpolate_marker_file);

//   // FILE *interpolate_oversampled_timestamps_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/interpolate_oversampled_float_timestamps.csv", "r");
//   // FILE *interpolate_values_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/interpolate_values.csv", "r");
//   // for (int i = 0; i < n; i++)
//   // {
//   //   got1 = fscanf(interpolate_oversampled_timestamps_file, "%lf", &interpolate_oversampled_float_timestamps[i]);
//   //   got2 = fscanf(interpolate_values_file, "%lf", &interpolate_values[i]);
//   //   if ((got1 + got2) != 2)
//   //     break; // wrong number of tokens - maybe end of file
//   //   if (interpolate_values[i] != 0)
//   //     continue;
//   // }
//   // fclose(interpolate_oversampled_timestamps_file);
//   // fclose(interpolate_values_file);

//   // FILE *fit_weighted_over_t_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/fit_weighted_over_t.csv", "r");
//   // FILE *fit_weighted_values_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/fit_weighted_values.csv", "r");
//   // FILE *fit_weighted_marker_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/fit_weighted_marker.csv", "r");
//   // for (int i = 0; i < n; i++)
//   // {
//   //   got1 = fscanf(fit_weighted_over_t_file, "%lf", &fit_weighted_over_t[i]);
//   //   got2 = fscanf(fit_weighted_values_file, "%lf", &fit_weighted_values[i]);
//   //   got3 = fscanf(fit_weighted_marker_file, "%lf", &fit_weighted_marker[i]);

//   //   if ((got1 + got2 + got3) != 3)
//   //     break; // wrong number of tokens - maybe end of file
//   //   if (fit_weighted_over_t[i] != 0)
//   //     continue;
//   // }
//   // fclose(fit_weighted_over_t_file);
//   // fclose(fit_weighted_values_file);
//   // fclose(fit_weighted_marker_file);

//   // FILE *mark_gamma_t_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/marker_gamma_t.csv", "r");
//   // FILE *mark_gamma_values_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/marker_gamma_values.csv", "r");
//   // for (int i = 0; i < n; i++)
//   // {
//   //   got1 = fscanf(mark_gamma_t_file, "%lf", &mark_gamma_t[i]);
//   //   got2 = fscanf(mark_gamma_values_file, "%lf", &mark_gamma_values[i]);
//   //   if ((got1 + got2) != 2)
//   //     break; // wrong number of tokens - maybe end of file
//   //   if (mark_gamma_t[i] != 0)
//   //     continue;
//   // }
//   // fclose(mark_gamma_t_file);
//   // fclose(mark_gamma_values_file);

//   FILE *single_col_t_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/single_col_t.csv", "r");
//   FILE *single_col_values_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/single_col_values.csv", "r");
//   for (int i = 0; i < n2; i++)
//   {
//     got1 = fscanf(single_col_t_file, "%lf", &single_col_t[i]);
//     got2 = fscanf(single_col_values_file, "%lf", &single_col_values[i]);
//     if ((got1 + got2) != 2)
//       break; // wrong number of tokens - maybe end of file
//     if (single_col_t[i] != 0)
//       continue;
//   }
//   fclose(single_col_t_file);
//   fclose(single_col_values_file);

//   // edges_t edges = extrapolate_edges(n, edges_marker, confident, sampling_rate);
//   // edges_t edges = extrapolate_neighbor(n, edges_marker, sampling_rate, k, confident);

//   // values_dataframe_t points_ex;
//   // points_ex.size = 1;
//   // points_ex.timestamps = malloc(sizeof(double));
//   // points_ex.timestamps[0] = 1485471758.48;
//   // points_ex.values = malloc(sizeof(double));
//   // points_ex.values[0] = -4.678013541955833;
//   // values_dataframe_t df = extrapolate_interpolate(n, interpolate_oversampled_float_timestamps, interpolate_values,
//   //                                                 interpolate_marker, points_ex.size, points_ex,
//   //                                                 sampling_rate, confident);

//   //  fit_weighted tests

//   //  left_start left_end  right_start right_end
//   //  15841      15845     15855       15851

//   //  fitted_left
//   //  nk min           range
//   //  7  1485471758.51 0.039999961853027344
//   //  coef
//   //  array([-3.9578835 , -3.87975805, -3.72350755, -3.48829143, -3.25230215, -3.09478554, -3.01602728])
//   //  knot
//   //  array([0.  , 0.  , 0.  , 0.  , 0.25, 0.5 , 0.75, 1.  , 1.  , 1.  , 1.  ])

//   //  fitted_right
//   //  nk min           range
//   //  7  1485471758.41 0.039999961853027344
//   //  coef
//   //  [-2.723240452872044, -2.8147348209162466, -2.9977221149619, -3.2754873149391814, -3.557491242420053, -3.7470792226270766, -3.841872750998506]
//   //  knot (11)
//   //  [0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0]
//   // smooth_spline_model_t fitted_model = fit_weighted(n, fit_weighted_over_t, fit_weighted_values,
//   //                                                   fit_weighted_marker, 15841, 15845, spar,
//   //                                                   sampling_rate, k);

//   // extrapolate_fit tests
//   // edges_t neighbors = extrapolate_neighbor(n, fit_weighted_marker, sampling_rate, k, confident);
//   // values_dataframe_t points_df = extrapolate_fit(n, fit_weighted_over_t, fit_weighted_values, neighbors,
//   //                                                fit_weighted_marker, spar, sampling_rate, k);

//   // mark_gamma tests
//   // noise 0.03
//   // range = -4, 4
//   // double *marker = mark_gamma(n, mark_gamma_t, mark_gamma_values, -4, 4, 0.03);

//   values_dataframe_t single_col_output = extrapolate_single_col(n2, single_col_t, single_col_values,
//                                                                 -4.0, 4.0, 0.03, 0.05, 0.6);

//   // FILE *single_result_values_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/single_col_output_values.csv", "a");
//   // FILE *single_result_timestamps_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/extrapolate/single_col_output_t.csv", "a");
//   // for (int i = 0; i < single_col_output.size; i++)
//   // {
//   //   fprintf(single_result_timestamps_file, "%f\n", single_col_output.timestamps[i]);
//   //   fprintf(single_result_values_file, "%f\n", single_col_output.values[i]);
//   // }
//   // fclose(single_result_timestamps_file);
//   // fclose(single_result_values_file);

//   run_subtests();
//   return 0;
// }
