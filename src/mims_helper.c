#include "mims_helper.h"

int get_break_size_in_seconds(int break_size, time_unit_t time_unit)
{
    int factor;
    switch (time_unit)
    {
    case second:
        factor = 1;
        break;
    case minute:
        factor = 60;
        break;
    case hour:
        factor = 3600;
        break;
    case day:
        factor = 86400; // 3600 * 24;
        break;
    }

    int break_size_in_seconds = break_size * factor;
    return break_size_in_seconds;
}

int parse_epoch_string(int break_size, time_unit_t time_unit, int sampling_rate)
{
    int num_seconds_in_epoch = (int)get_break_size_in_seconds(break_size, time_unit);
    int num_samples = num_seconds_in_epoch * sampling_rate;
    return num_samples;
}

int get_sampling_rate(dataframe_t *dataframe)
{
    int duration_in_seconds = (int)(dataframe->timestamps[dataframe->size - 1] - dataframe->timestamps[0]);
    int sampling_rate = dataframe->size / duration_in_seconds;
    return sampling_rate;
}

static void compute_time_segments(
    dataframe_t *dataframe, double start_time, double break_size_in_seconds)
{
    dataframe->n_segments = (int)((dataframe->timestamps[dataframe->size - 1] - start_time) / break_size_in_seconds) + 1;
    dataframe->segments = malloc(dataframe->n_segments * sizeof(int));

    int segment_i = 0;
    int last_break_index = 0;
    for (int i = 0; i < dataframe->size; i++)
    {
        if (((dataframe->timestamps[i] - dataframe->timestamps[last_break_index]) >= break_size_in_seconds) || (i == dataframe->size - 1))
        {
            dataframe->segments[segment_i] = last_break_index;
            last_break_index = i;
            segment_i += 1;
        }
    }
    return;
}

// Segments the input sensor dataframe into epoch windows with length specified in break_size.
// Mutates
void segment_data(dataframe_t *dataframe, int break_size, time_unit_t time_unit, double start_time)
{
    double break_size_in_seconds = (double)get_break_size_in_seconds(break_size, time_unit);

    compute_time_segments(dataframe, start_time, break_size_in_seconds);

    return;
}

int sequence_length(double start, double stop, double step)
{
    return (int)(((stop - start) / step) + pow(1, -10));
}

// double *sequence(double start, double stop, double step)
// {
//     int sequence_len = sequence_length(start, stop, step);
//     double *sequence = malloc(sequence_len * sizeof(double));

//     double current = start;
//     for (int i = 0; i < sequence_len; i++)
//     {
//         sequence[i] = current;
//         current += step;
//     }

//     return sequence;
// }

double *sequence(double start, double stop, double step)
{
    int sequence_len = sequence_length(start, stop, step);
    double *sequence = malloc(sequence_len * sizeof(double));

    sequence[0] = start;
    // sequence[sequence_len - 1] = stop;
    for (int i = 1; i < sequence_len; i++)
    {
        sequence[i] = start + (double)(i * step);
    }
    return sequence;
}

double *linspace(double start, double stop, int n)
{
    double *sequence = malloc(n * sizeof(double));
    double step = (stop - start) / (double)(n - 1.0);
    double current = start;
    for (int i = 0; i < n; i++)
    {
        sequence[i] = current;
        current += step;
    }
    return sequence;
}
