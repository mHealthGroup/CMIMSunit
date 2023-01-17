#include "basic_tests.h"
#include "../src/mims_unit.h"

int main(int argc, char **argv)
{
    // Data Set 1
    char input_filename[66] = "/Users/arytonhoi/Kode/mhealth/cmims/data/mims_unit/test_2/raw.csv";
    char expected_output_filename[71] = "/Users/arytonhoi/Kode/mhealth/cmims/data/mims_unit/test_2/r_output.csv";

    // Data Set 2
    // char input_filename[72] = "/Users/arytonhoi/Kode/mhealth/cmims/data/mims_unit/aditya_sleep/raw.csv";
    // char expected_output_filename[77] = "/Users/arytonhoi/Kode/mhealth/cmims/data/mims_unit/aditya_sleep/r_output.csv";

    dataframe_t *mims_data = consistency_test(input_filename);
    precision_test(mims_data, expected_output_filename);
    before_after_df_test(input_filename, expected_output_filename);
    measure_runtime(input_filename);

    free_dataframe(mims_data);

    return 0;
}
