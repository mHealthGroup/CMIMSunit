#include <time.h>
#include "basic_tests.h"

dataframe_t consistency_test(char *input_filename)
{
    dataframe_t input_df = read_csv(input_filename);
    dataframe_t previous_output, current_output;
    previous_output = mims_unit(&input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1);
    for (int i = 0; i < 20; i++)
    {
        current_output = (i % 2)
                             ? mims_unit(&input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1)
                             : mims_unit_from_filename(input_filename, -8, 8, 1, minute, 0.03, 0.05,
                                                       0.6, 0.2, 5.0, 1);

        for (int j = 0; j < previous_output.size; j++)
        {
            if (previous_output.mims_data[j] != current_output.mims_data[j])
            {
                printf("Failed consistency_test\n");
                return current_output;
            }
        }
        previous_output = current_output;
    }

    printf("Passed consistency_test\n");
    return current_output;
}

void precision_test(dataframe_t *output_df, char *expected_output_filename)
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
            printf("Failed precision test\n");
            return;
        }
    }

    printf("Passed precision test\n");
    return;
}

void before_after_df_test(char *input_filename, char *expected_output_filename)
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

    dataframe_t mims_data = custom_mims_unit(&middle_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2,
                                             5.0, 1, &before_df, &after_df);
    precision_test(&mims_data, expected_output_filename);
}

void measure_runtime(char *input_filename)
{
    dataframe_t input_df = read_csv(input_filename);

    clock_t begin = clock();
    mims_unit(&input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1);
    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Cmims ran in %f seconds\n", time_spent);
    return;
}
