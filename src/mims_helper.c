#include "mims_helper.h"

uint32_t get_break_size_in_seconds(uint16_t break_size, time_unit_t time_unit)
{
    uint32_t factor;
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

    uint32_t break_size_in_seconds = break_size * factor;
    return break_size_in_seconds;
}

uint32_t parse_epoch_string(uint16_t break_size, time_unit_t time_unit, uint16_t sampling_rate)
{
    uint32_t num_seconds_in_epoch = get_break_size_in_seconds(break_size, time_unit);
    uint32_t num_samples = num_seconds_in_epoch * sampling_rate;
    return num_samples;
}

uint16_t get_sampling_rate(dataframe_t *dataframe)
{
    uint32_t duration_in_seconds = (uint32_t)(dataframe->timestamps[dataframe->size - 1] - dataframe->timestamps[0]);
    uint16_t sampling_rate = dataframe->size / duration_in_seconds;
    return sampling_rate;
}

static void compute_time_segments(
    dataframe_t *dataframe, double start_time, uint32_t break_size_in_seconds)
{
    dataframe->n_segments = (uint32_t)((dataframe->timestamps[dataframe->size - 1] - start_time) / (double)break_size_in_seconds) + 1;
    dataframe->segments = malloc(dataframe->n_segments * sizeof(uint32_t));

    uint32_t segment_i = 0;
    uint32_t last_break_index = 0;
    for (uint32_t i = 0; i < dataframe->size; i++)
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
void segment_data(dataframe_t *dataframe, uint16_t break_size, time_unit_t time_unit, double start_time)
{
    uint32_t break_size_in_seconds = get_break_size_in_seconds(break_size, time_unit);
    compute_time_segments(dataframe, start_time, break_size_in_seconds);
    return;
}

uint32_t sequence_length(double start, double stop, double step)
{
    return (uint32_t)(((stop - start) / step) + pow(1, -10));
}

double *sequence(double start, double stop, double step)
{
    uint32_t sequence_len = sequence_length(start, stop, step);
    double *sequence = malloc(sequence_len * sizeof(double));

    sequence[0] = start;
    for (uint32_t i = 1; i < sequence_len; i++)
    {
        sequence[i] = start + (double)(i * step);
    }

    return sequence;
}

double *linspace(double start, double stop, uint32_t n)
{
    double *sequence = malloc(n * sizeof(double));
    double step = (stop - start) / (double)(n - 1.0);
    double current = start;
    for (uint32_t i = 0; i < n; i++)
    {
        sequence[i] = current;
        current += step;
    }

    return sequence;
}
