#include "config.h"
#include <stdio.h>
#include "mims_helper.h"
#include "mims_unit.h"

void show_version()
{
    printf("CMIMS Version: %d.%d\n", CMIMS_VERSION_MAJOR, CMIMS_VERSION_MINOR);
}

void print_dataframe(dataframe_t *df, int head)
{
    for (int i = 0; i < head; i++)
    {
        printf("%f,%f,%f,%f\n", df->x[i], df->y[i], df->z[i], df->mims_data[i]);
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        // report version
        show_version();
    }
    else
    {
        char *filepath = argv[1];

        dataframe_t input_df = read_csv(filepath);
        dataframe_t output_df = mims_unit(&input_df, -8, 8, 1, minute, 0.03, 0.05, 0.6, 0.2, 5.0, 1);
        print_dataframe(&output_df, 10);
    }
    return 0;
}