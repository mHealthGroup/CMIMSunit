#include <time.h>
#include "basic_tests.h"
// #include "../src/mims_helper.h"

dataframe_t *consistency_test(char *input_filename)
{
    dataframe_t *input_df = read_csv(input_filename);
    dataframe_t *previous_output, *current_output;
    previous_output = mims_unit(input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1);
    for (int i = 0; i < 20; i++)
    {
        current_output = (i % 2)
                             ? mims_unit(input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1)
                             : mims_unit_from_filename(input_filename, -8, 8, 1, minute, 0.03, 0.05,
                                                       0.6, 0.2, 5.0, 1);

        for (int j = 0; j < previous_output->size; j++)
        {
            if (previous_output->mims_data[j] != current_output->mims_data[j])
            {
                printf("Failed consistency_test\n");
                return current_output;
            }
        }
        free_dataframe(previous_output);
        previous_output = current_output;
    }

    printf("Passed consistency_test\n");
    free_dataframe(input_df);
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
            free(r_output);
            return;
        }
    }

    printf("Passed precision test\n");
    free(r_output);
    return;
}

void before_after_df_test(char *input_filename, char *expected_output_filename)
{
    dataframe_t *input_df = read_csv(input_filename);

    dataframe_t *before_df = create_dataframe(
        1,
        malloc(sizeof(double)),
        malloc(sizeof(double)),
        malloc(sizeof(double)),
        malloc(sizeof(double)),
        0, NULL, NULL);
    dataframe_t *after_df = create_dataframe(
        1,
        malloc(sizeof(double)),
        malloc(sizeof(double)),
        malloc(sizeof(double)),
        malloc(sizeof(double)),
        0, NULL, NULL);
    before_df->timestamps[0] = input_df->timestamps[0];
    before_df->x[0] = input_df->x[0];
    before_df->y[0] = input_df->y[0];
    before_df->z[0] = input_df->z[0];

    after_df->timestamps[0] = input_df->timestamps[input_df->size - 1];
    after_df->x[0] = input_df->x[input_df->size - 1];
    after_df->y[0] = input_df->y[input_df->size - 1];
    after_df->z[0] = input_df->z[input_df->size - 1];

    dataframe_t *middle_df = create_dataframe(
        input_df->size - 2,
        malloc((input_df->size - 2) * sizeof(double)),
        malloc((input_df->size - 2) * sizeof(double)),
        malloc((input_df->size - 2) * sizeof(double)),
        malloc((input_df->size - 2) * sizeof(double)),
        0, NULL, NULL);

    for (int i = 0; i < middle_df->size; i++)
    {
        middle_df->timestamps[i] = input_df->timestamps[i + 1];
        middle_df->x[i] = input_df->x[i + 1];
        middle_df->y[i] = input_df->y[i + 1];
        middle_df->z[i] = input_df->z[i + 1];
    }

    dataframe_t *mims_data = custom_mims_unit_before_after_dataframe(middle_df, -8, 8, 1, minute,
                                                                     0.03, 0.05, 0.6, 0.2, 5.0, 1,
                                                                     before_df, after_df);
    precision_test(mims_data, expected_output_filename);
    free_dataframe(input_df);
    free_dataframe(before_df);
    free_dataframe(after_df);
    free_dataframe(middle_df);
    free_dataframe(mims_data);
}

void measure_runtime(char *input_filename)
{
    dataframe_t *input_df = read_csv(input_filename);

    clock_t begin = clock();
    dataframe_t *output_df = mims_unit(input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1);
    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Cmims ran in %f seconds\n", time_spent);

    free_dataframe(input_df);
    free_dataframe(output_df);
    return;
}

void test_sm_spline()
{
    int n = 5;
    double spar = 0.6;
    double over_t[5] = {
        1485471758.410000085830688476562500000000000000000000000000000000000000000000000000000000000,
        1485471758.420000076293945312500000000000000000000000000000000000000000000000000000000000000,
        1485471758.430000066757202148437500000000000000000000000000000000000000000000000000000000000,
        1485471758.440000057220458984375000000000000000000000000000000000000000000000000000000000000,
        1485471758.450000047683715820312500000000000000000000000000000000000000000000000000000000000,
    };
    double over_values[5] = {
        -2.72324045287203775345119538542348891496658325195312500000000000000000000000000000000000000,
        -2.99826949144669985258815358974970877170562744140625000000000000000000000000000000000000000,
        -3.27619376952312002515554922865703701972961425781250000000000000000000000000000000000000000,
        -3.55788758289167850179524066334124654531478881835937500000000000000000000000000000000000000,
        -3.84187275099852998394567293871659785509109497070312500000000000000000000000000000000000000,
    };
    double weights[5] = {
        0.849863384924921416718746058904798701405525207519531250000000000000000000000000000000000000,
        0.849863384924921416718746058904798701405525207519531250000000000000000000000000000000000000,
        0.849863384924921416718746058904798701405525207519531250000000000000000000000000000000000000,
        0.849863384924921416718746058904798701405525207519531250000000000000000000000000000000000000,
        1.600546460300313889035805914318189024925231933593750000000000000000000000000000000000000000,
    };

    double expected[5] = {
        -2.720145075829855763061,
        -2.999651194341386162279,
        -3.279436632826645769967,
        -3.559638255769206338641,
        -3.840131199147975848973};

    smooth_spline_model_t *model = sm_spline_coef(n, over_t, over_values, n, weights, spar);
    double *results = predict_smooth_spline(model, over_t, n, 0);
    free_smooth_spline_model(model);

    for (int i = 0; i < n; i++)
    {
        if (fabs(results[i] - expected[i]) > 0)
        {
            printf("Failed test_sm_spline\n");
            return;
        }
    }

    free(results);
    printf("Pased test_sm_spline\n");
}