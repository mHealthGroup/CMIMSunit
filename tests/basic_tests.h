#ifndef _BASIC_TESTS_H_
#define _BASIC_TESTS_H_

#include "../src/mims_unit.h"

dataframe_t consistency_test(char *input_filename);
void precision_test(dataframe_t *output_df, char *expected_output_filename);
void before_after_df_test(char *input_filename, char *expected_output_filename);
void measure_runtime(char *input_filename);

#endif // _BASIC_TESTS_H_
