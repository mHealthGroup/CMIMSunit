#ifndef _STATS_SMSPLINE_H_
#define _STATS_SMSPLINE_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct
{
    double *coef;
    double *ty;
    double *lev;
    double spar;
    double *parms;
    double *crit;
    int *iparms;
    int *ier;
    double *scrtch;
} fit_t;

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
    uint8_t no_weights;
    double *lev;
    double cv_crit;
    // double pen_crit; // Unused
    double *crit;
    double df;
    double spar;
    double ratio;
    double lambda;
    // int *iparms; // Unused
    fit_t *fit;
    double *fit_knot;
    int fit_nk;
    double fit_min;
    double fit_range;
    double *fit_coef;
} smooth_spline_model_t;

smooth_spline_model_t *sm_spline_coef(const int n, double *x, double *y, const int w_len, double *w,
                                      double spar);
double *predict_smooth_spline(const smooth_spline_model_t *model, const double *x, const int x_len,
                              const int deriv);
void free_smooth_spline_model(smooth_spline_model_t *model);

#endif // _STATS_SMSPLINE_H_