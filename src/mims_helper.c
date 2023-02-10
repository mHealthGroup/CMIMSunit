#include "mims_helper.h"
#include <string.h> // strtok

#define LINE_COUNT_BUFFER_SIZE 65535
#define CSV_LINE_BUFFER_SIZE 128

static uint32_t get_break_size_in_seconds(const uint16_t break_size, const time_unit_t time_unit)
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

static void compute_time_segments(
    dataframe_t *dataframe, const double start_time, const uint32_t break_size_in_seconds)
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

uint32_t parse_epoch_string(const uint16_t break_size, const time_unit_t time_unit,
                            const uint16_t sampling_rate)
{
    uint32_t num_seconds_in_epoch = get_break_size_in_seconds(break_size, time_unit);
    uint32_t num_samples = num_seconds_in_epoch * sampling_rate;
    return num_samples;
}

uint16_t get_sampling_rate(const dataframe_t *dataframe)
{
    uint32_t duration_in_seconds = (uint32_t)(dataframe->timestamps[dataframe->size - 1] - dataframe->timestamps[0]);
    uint16_t sampling_rate = dataframe->size / duration_in_seconds;
    return sampling_rate;
}

// Segments the input sensor dataframe into epoch windows with length specified in break_size.
void segment_data(dataframe_t *dataframe, const uint16_t break_size,
                  const time_unit_t time_unit, const double start_time)
{
    uint32_t break_size_in_seconds = get_break_size_in_seconds(break_size, time_unit);
    compute_time_segments(dataframe, start_time, break_size_in_seconds);
    return;
}

uint32_t sequence_length(const double start, const double stop, const double step)
{
    return (uint32_t)(((stop - start) / step) + pow(1, -10));
}

double *sequence(const double start, const double stop, const double step)
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

double *linspace(const double start, const double stop, const uint32_t n)
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

int count_lines(const char *filename)
{
    FILE *file = fopen(filename, "r");
    char buf[LINE_COUNT_BUFFER_SIZE];
    int counter = 0;
    int i;
    for (;;)
    {
        size_t res = fread(buf, 1, LINE_COUNT_BUFFER_SIZE, file);
        if (ferror(file))
        {
            fclose(file);
            return -1;
        }

        for (i = 0; i < res; i++)
            if (buf[i] == '\n')
                counter++;

        if (feof(file))
            break;
    }

    fclose(file);
    return counter;
}

dataframe_t *read_csv(const char *filename)
{
    int n = count_lines(filename) - 1; // Exclude header row
    FILE *file = fopen(filename, "r");
    // rewind(file);
    dataframe_t *df = create_dataframe(
        n,
        malloc(n * sizeof(double)),
        malloc(n * sizeof(double)),
        malloc(n * sizeof(double)),
        malloc(n * sizeof(double)),
        0, NULL, NULL);

    char line_buffer[CSV_LINE_BUFFER_SIZE];
    char *token;
    for (int i = -1; i < n; i++)
    {
        if (!fgets(line_buffer, CSV_LINE_BUFFER_SIZE, file))
            break;
        if (i == -1)
            continue; // skip first line

        token = strtok(line_buffer, ",");
        df->timestamps[i] = atof(token);
        token = strtok(NULL, ",");
        df->x[i] = atof(token);
        token = strtok(NULL, ",");
        df->y[i] = atof(token);
        token = strtok(NULL, ",");
        df->z[i] = atof(token);
    }

    fclose(file);
    return df;
}

dataframe_t *create_dataframe(uint32_t size, double *timestamps, double *x,
                              double *y, double *z, uint32_t n_segments,
                              uint32_t *segments, double *mims_data)
{
    dataframe_t *dataframe = malloc(sizeof(dataframe_t));

    dataframe->size = size;
    dataframe->timestamps = timestamps;
    dataframe->x = x;
    dataframe->y = y;
    dataframe->z = z;
    dataframe->n_segments = n_segments;
    dataframe->segments = segments;
    dataframe->mims_data = mims_data;

    return dataframe;
}

void free_dataframe(dataframe_t *df)
{
    free(df->timestamps);
    free(df->x);
    free(df->y);
    free(df->z);
    free(df->segments);
    free(df->mims_data);
    free(df);

    return;
}

dataframe_t *concat_dataframes(const dataframe_t *df_1, const dataframe_t *df_2)
{
    uint32_t n = df_1->size + df_2->size;
    double *timestamps = malloc(n * sizeof(double));
    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));
    double *z = malloc(n * sizeof(double));

    memcpy(timestamps, df_1->timestamps, df_1->size * sizeof(*df_1->timestamps));
    memcpy(x, df_1->x, df_1->size * sizeof(*df_1->x));
    memcpy(y, df_1->y, df_1->size * sizeof(*df_1->y));
    memcpy(z, df_1->z, df_1->size * sizeof(*df_1->z));

    memcpy(timestamps + df_1->size, df_2->timestamps, df_2->size * sizeof(*df_2->timestamps));
    memcpy(x + df_1->size, df_2->x, df_2->size * sizeof(*df_2->x));
    memcpy(y + df_1->size, df_2->y, df_2->size * sizeof(*df_2->y));
    memcpy(z + df_1->size, df_2->z, df_2->size * sizeof(*df_2->z));

    dataframe_t *output_df = create_dataframe(n, timestamps, x, y, z, 0, NULL, NULL);
    return output_df;
}
