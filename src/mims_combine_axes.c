#include "mims_combine_axes.h"

double *sum_up(dataframe_t *dataframe)
{
  int n = dataframe->size;
  double *output = malloc(n * sizeof(double));
  for (int i = 0; i < n; i++)
    output[i] = dataframe->x[i] + dataframe->y[i] + dataframe->z[i];

  return output;
}

// int main(int argc, char **argv)
// {
//   int n = 108000;
//   int m = 720;

//   int got1, got2, got3, got4;

//   double *timestamps = malloc(n * sizeof(double));
//   double *x = malloc(n * sizeof(double));
//   double *y = malloc(n * sizeof(double));
//   double *z = malloc(n * sizeof(double));
//   int *segments = malloc(m * sizeof(int));
//   FILE *t_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/aggregate/integrate_mims_timestamps.csv", "r");
//   FILE *x_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/aggregate/integrate_mims_x.csv", "r");
//   FILE *y_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/aggregate/integrate_mims_y.csv", "r");
//   FILE *z_file = fopen("/Users/arytonhoi/Kode/mhealth/cmims/test/data/aggregate/integrate_mims_z.csv", "r");
//   for (int i = 0; i < n; i++)
//   {
//     got1 = fscanf(t_file, "%lf", &timestamps[i]);
//     timestamps[i] = timestamps[i] / pow(10, 9);
//     got2 = fscanf(x_file, "%lf", &x[i]);
//     got3 = fscanf(y_file, "%lf", &y[i]);
//     got4 = fscanf(z_file, "%lf", &z[i]);
//     if ((got1 + got2 + got3 + got4) != 4)
//       break; // wrong number of tokens - maybe end of file
//     if (x[i] != 0)
//       continue;
//   }
//   fclose(t_file);
//   fclose(x_file);
//   fclose(y_file);
//   fclose(z_file);

//   dataframe_t input_df;
//   input_df.size = n;
//   input_df.timestamps = timestamps;
//   input_df.x = x;
//   input_df.y = y;
//   input_df.z = z;
//   double *result = sum_up(&input_df);

//   return 0;
// }
