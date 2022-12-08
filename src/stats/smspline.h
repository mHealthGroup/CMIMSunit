#ifndef _STATS_SMSPLINE_H_
#define _STATS_SMSPLINE_H_

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct
{
    double *x;
    double *y;
    double *w;
    double *yin;
    double tol;
    double *data_x;
    double *data_y;
    double *data_w;
    bool no_weights;
    double *lev;
    double cv_crit;
    // double pen_crit;
    double *crit;
    double df;
    double spar;
    double ratio;
    double lambda;
    // int *iparms;
    double *fit_knot;
    int fit_nk;
    double fit_min;
    double fit_range;
    double *fit_coef;
} smooth_spline_model_t;

smooth_spline_model_t SmSplineCoef(int n, double *x, double *y, int w_len, double *w, double spar);
double *predict_smooth_spline(smooth_spline_model_t model, double *x, int x_len, int deriv);

#endif // _STATS_SMSPLINE_H_