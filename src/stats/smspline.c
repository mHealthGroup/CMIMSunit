#include "smspline.h"

typedef struct
{
    double low;
    double high;
    double tol;
    double eps;
    int maxit;
    uint8_t trace;
} contr_sp_t;

static double iqr(int n, double *x)
{
    float index[2] = {1 + (n - 1) * 0.25, 1 + (n - 1) * 0.75};
    int lo[2] = {(int)floor(index[0]), (int)floor(index[1])};
    int hi[2] = {(int)ceil(index[0]), (int)ceil(index[1])};

    double qs[2] = {x[lo[0]], x[lo[1]]};

    // R code, probably unreached
    // i <- which(index > lo)
    // h <- (index - lo)[i] # > 0	by construction
    // qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]

    return qs[1] - qs[0];
}

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

static int _nknots_smspl(int n)
{
    if (n < 50)
    {
        return n;
    }
    else
    {
        double a1 = log2(50);
        double a2 = log2(100);
        double a3 = log2(140);
        double a4 = log2(200);
        if (n < 200)
            return (int)pow(2, (a1 + (a2 - a1) * (n - 50) / 150));
        else if (n < 800)
            return (int)pow(2, (a2 + (a3 - a2) * (n - 200) / 600));
        else if (n < 3200)
            return (int)pow(2, (a3 + (a4 - a3) * (n - 800) / 2400));
        else
            return 200 + pow(n - 3200, 0.2);
    }
}

static double weighted_avg(int n, double *x, double *weights)
{
    double weight_sum = 0;
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += x[i] * weights[i];
        weight_sum += weights[i];
    }

    return sum / weight_sum;
}

static double mean(int n, double *x)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += x[i];
    }
    return sum / (double)n;
}

static double fsign(double x, double y)
{
    return ((y >= 0) ? fabs(x) : -fabs(x));
}

static int convert_2d_to_1d_idx(int num_rows, int row_num, int col_num)
{
    row_num = row_num + 1;
    col_num = col_num + 1;
    return ((col_num - 1) * num_rows + row_num) - 1;
}

static int findInterval2(double *xt, int n, double x,
                         uint8_t rightmost_closed, uint8_t all_inside,
                         uint8_t left_open, // <- new in findInterval2()
                         int ilo, int *mflag)
{
    int istep, middle, ihi;

#define left_boundary                                                      \
    {                                                                      \
        *mflag = -1;                                                       \
        return ((all_inside || (rightmost_closed && x == xt[1])) ? 1 : 0); \
    }

#define right_boundary                                           \
    {                                                            \
        *mflag = +1;                                             \
        return ((all_inside || (rightmost_closed && x == xt[n])) \
                    ? (n - 1)                                    \
                    : n);                                        \
    }

#define X_grtr(XT_v) x > XT_v || (!left_open && x >= XT_v)
#define X_smlr(XT_v) x < XT_v || (left_open && x <= XT_v)

    if (n == 0)
    {
        *mflag = 0;
        return 0;
    }
    /* 1-indexing : */
    --xt;

    if (ilo <= 0)
    {
        if (X_smlr(xt[1]))
            left_boundary;
        ilo = 1;
    }
    ihi = ilo + 1;
    if (ihi >= n)
    {
        if (X_grtr(xt[n]))
            right_boundary;
        if (n <= 1) /* x < xt[1] */
            left_boundary;
        ilo = n - 1;
        ihi = n;
    }

    if (X_smlr(xt[ihi]))
    {
        if (X_grtr(xt[ilo]))
        {
            /* `lucky': same interval as last time */
            *mflag = 0;
            return ilo;
        }
        /* **** now x < xt[ilo] .	decrease  ilo  to capture  x */
        if (!left_open)
            for (istep = 1;; istep *= 2)
            {
                ihi = ilo;
                ilo = ihi - istep;
                if (ilo <= 1)
                    break;
                if (x >= xt[ilo])
                    goto L50;
            }
        else
            for (istep = 1;; istep *= 2)
            {
                ihi = ilo;
                ilo = ihi - istep;
                if (ilo <= 1)
                    break;
                if (x > xt[ilo])
                    goto L51;
            }
        ilo = 1;
        if (X_smlr(xt[1]))
            left_boundary;
    }
    else
    {
        /* **** now x >= xt[ihi] .	increase  ihi  to capture  x */
        if (!left_open)
            for (istep = 1;; istep *= 2)
            {
                ilo = ihi;
                ihi = ilo + istep;
                if (ihi >= n)
                    break;
                if (x < xt[ihi])
                    goto L50;
            }
        else
            for (istep = 1;; istep *= 2)
            {
                ilo = ihi;
                ihi = ilo + istep;
                if (ihi >= n)
                    break;
                if (x <= xt[ihi])
                    goto L51;
            }
        if (X_grtr(xt[n]))
            right_boundary;
        ihi = n;
    }

    if (left_open)
        goto L51; /* There _is_ a path to here, avoiding return and goto */

L50: // ! left_open
    /* **** now xt[ilo] <= x < xt[ihi] . narrow the interval. */
    for (;;)
    {
        middle = (ilo + ihi) / 2;
        if (middle == ilo)
        {
            *mflag = 0;
            return ilo;
        }
        /* note. it is assumed that middle = ilo in case ihi = ilo+1 . */
        if (x >= xt[middle])
            ilo = middle;
        else
            ihi = middle;
    }

L51: // left_open
    /* **** now xt[ilo] < x <= xt[ihi] . narrow the interval. */
    for (;;)
    {
        middle = (ilo + ihi) / 2;
        if (middle == ilo)
        {
            *mflag = 0;
            return ilo;
        }
        /* note. it is assumed that middle = ilo in case ihi = ilo+1 . */
        if (x > xt[middle])
            ilo = middle;
        else
            ihi = middle;
    }
} /* findInterval2 */

// has been in API -- keep for compatibility:
static int findInterval(double *xt, int n, double x,
                        uint8_t rightmost_closed, uint8_t all_inside, int ilo,
                        int *mflag)
{
    return findInterval2(xt, n, x, rightmost_closed, all_inside, 0, ilo, mflag);
}

static int interv(double *xt, int n, double x,
                  uint8_t rightmost_closed, uint8_t all_inside,
                  int ilo, int *mflag)
{
    return findInterval(xt, n, x, rightmost_closed, all_inside, ilo, mflag);
}

static void bsplvb(double *t, int lent, int jhigh, int index, double x, int left, double *biatx,
                   int *j, double *deltal, double *deltar)
{
    int i;

    if (index != 2)
    {
        *j = 1;
        biatx[1 - 1] = 1.0;
        if (*j >= jhigh)
            return;
    }

    while (1)
    {
        int jp1 = *j + 1;
        deltar[*j - 1] = t[left + *j - 1] - x;
        deltal[*j - 1] = x - t[left + 1 - *j - 1];
        double saved = 0.0;
        for (i = 1; i < *j + 1; i++)
        {
            double term = biatx[i - 1] / (deltar[i - 1] + deltal[jp1 - i - 1]);
            biatx[i - 1] = saved + deltar[i - 1] * term;
            saved = deltal[jp1 - i - 1] * term;
        }
        biatx[jp1 - 1] = saved;
        *j = jp1;
        if (!(*j < jhigh))
            break;
    }
}

static void bsplvd(double *t, int lent, int k, double x, int left, double *a, double *dbiatx, int nderiv,
                   int *j, double *deltal, double *deltar)
{
    int mhigh = MAX(MIN(nderiv, k), 1);
    int kp1 = k + 1;

    bsplvb(t, lent, kp1 - mhigh, 1, x, left, dbiatx, j, deltal, deltar);

    if (mhigh == 1)
        return;

    int ideriv = mhigh;
    int m_idx, j_idx, i_idx, l_idx;
    for (m_idx = 2; m_idx < mhigh + 1; m_idx++)
    {
        int jp1mid = 1;
        for (j_idx = ideriv; j_idx < k + 1; j_idx++)
        {
            dbiatx[convert_2d_to_1d_idx(4, j_idx - 1, ideriv - 1)] = dbiatx[convert_2d_to_1d_idx(4, jp1mid - 1, 1 - 1)];
            jp1mid += 1;
        }

        ideriv -= 1;
        bsplvb(t, lent, kp1 - ideriv, 2, x, left, dbiatx, j, deltal, deltar);
    }

    int jlow = 1;
    for (i_idx = 1; i_idx < k + 1; i_idx++)
    {
        for (j_idx = jlow; j_idx < k + 1; j_idx++)
            a[convert_2d_to_1d_idx(4, j_idx - 1, i_idx - 1)] = 0.0;

        jlow = i_idx;
        a[convert_2d_to_1d_idx(4, i_idx - 1, i_idx - 1)] = 1.0;
    }

    for (m_idx = 2; m_idx < mhigh + 1; m_idx++)
    {
        int kp1mm = kp1 - m_idx;
        double fkp1mm = (double)kp1mm;
        int il = left;
        i_idx = k;

        for (l_idx = 1; l_idx < kp1mm + 1; l_idx++)
        {
            double factor = fkp1mm / (t[il + kp1mm - 1] - t[il - 1]);
            // the assumption that t(left) < t(left+1) makes denominator in  factor  nonzero.
            for (j_idx = 1; j_idx < i_idx + 1; j_idx++)
            {
                a[convert_2d_to_1d_idx(4, i_idx - 1, j_idx - 1)] = (a[convert_2d_to_1d_idx(4, i_idx - 1, j_idx - 1)] - a[convert_2d_to_1d_idx(4, i_idx - 1 - 1, j_idx - 1)]) * factor;
            }
            il -= 1;
            i_idx -= 1;
        }

        for (i_idx = 1; i_idx < k + 1; i_idx++)
        {
            double _sum = 0.0;
            int jlow = MAX(i_idx, m_idx);
            for (j_idx = jlow; j_idx < k + 1; j_idx++)
                _sum += a[convert_2d_to_1d_idx(4, j_idx - 1, i_idx - 1)] * dbiatx[convert_2d_to_1d_idx(4, j_idx - 1, m_idx - 1)];

            dbiatx[convert_2d_to_1d_idx(4, i_idx - 1, m_idx - 1)] = _sum;
        }
    }
}

static void sgram(double *sg0, double *sg1, double *sg2, double *sg3, double *tb, int nb,
                  int *j, double *deltal, double *deltar)
{
    int i, ii, jj;
    int mflag = 0;

    double *vnikx = malloc(12 * sizeof(double));
    double *work = malloc(16 * sizeof(double));
    double *yw1 = malloc(4 * sizeof(double));
    double *yw2 = malloc(4 * sizeof(double));

    int lentb = nb + 4;
    //  Initialise the sigma vectors
    for (i = 1; i < nb + 1; i++)
    {
        sg0[i - 1] = 0.0;
        sg1[i - 1] = 0.0;
        sg2[i - 1] = 0.0;
        sg3[i - 1] = 0.0;
    }

    int ileft = 1;
    for (i = 1; i < nb + 1; i++)
    {
        // Calculate a linear approximation to the second derivative of the
        // non-zero B-splines over the interval [tb(i),tb(i+1)].

        ileft = interv(
            (double *)tb,
            nb + 1,
            tb[i - 1],
            0,
            0,
            ileft,
            &mflag);

        // Left end second derivatives
        bsplvd(
            tb,
            lentb,
            4,
            tb[i - 1],
            ileft,
            work,
            vnikx,
            3,
            j,
            deltal,
            deltar);

        // Put values into yw1
        for (ii = 1; ii < 4 + 1; ii++)
        {
            yw1[ii - 1] = vnikx[convert_2d_to_1d_idx(4, ii - 1, 3 - 1)];
        }

        // Right end second derivatives
        bsplvd(
            tb, lentb, 4, tb[i + 1 - 1], ileft, work, vnikx, 3, j, deltal, deltar);

        // Slope*(length of interval) in Linear Approximation to B''
        for (ii = 1; ii < 4 + 1; ii++)
        {
            yw2[ii - 1] = vnikx[convert_2d_to_1d_idx(4, ii - 1, 3 - 1)] - yw1[ii - 1];
        }

        // Calculate Contributions to the sigma vectors
        double wpt = tb[i + 1 - 1] - tb[i - 1];
        if (ileft >= 4)
        {
            for (ii = 1; ii < 4 + 1; ii++)
            {
                jj = ii;
                sg0[ileft - 4 + ii - 1] = sg0[ileft - 4 + ii - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);

                jj = ii + 1;
                if (jj <= 4)
                {
                    sg1[ileft + ii - 4 - 1] = sg1[ileft + ii - 4 - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);
                }

                jj = ii + 2;
                if (jj <= 4)
                {
                    sg2[ileft + ii - 4 - 1] = sg2[ileft + ii - 4 - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);
                }

                jj = ii + 3;
                if (jj <= 4)
                {
                    sg3[ileft + ii - 4 - 1] = sg3[ileft + ii - 4 - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);
                }
            }
        }
        else if (ileft == 3)
        {
            for (ii = 1; ii < 3 + 1; ii++)
            {
                jj = ii;
                sg0[ileft - 3 + ii - 1] = sg0[ileft - 3 + ii - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);

                jj = ii + 1;
                if (jj <= 3)
                {
                    sg1[ileft + ii - 3 - 1] = sg1[ileft + ii - 3 - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);
                }

                jj = ii + 2;
                if (jj <= 3)
                {
                    sg2[ileft + ii - 3 - 1] = sg2[ileft + ii - 3 - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);
                }
            }
        }

        else if (ileft == 2)
        {
            for (ii = 1; ii < 2 + 1; i++)
            {
                jj = ii;
                sg0[ileft - 2 + ii - 1] = sg0[ileft - 2 + ii - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);

                jj = ii + 1;
                if (jj <= 2)
                {
                    sg1[ileft + ii - 2 - 1] = sg1[ileft + ii - 2 - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);
                }
            }
        }

        else if (ileft == 1)
        {
            ii = 1;
            jj = ii;
            sg0[ileft - 1 + ii - 1] = sg0[ileft - 1 + ii - 1] + wpt * (yw1[ii - 1] * yw1[jj - 1] + (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1]) * 0.5 + yw2[ii - 1] * yw2[jj - 1] * 0.3330);
        }
    }

    free(yw1);
    free(yw2);
    free(vnikx);
    free(work);
    return;
}

static void stxwx(
    double *x, double *z, double *w, int k, double *xknot, int n,
    double *y, double *hs0, double *hs1, double *hs2, double *hs3,
    int *j, double *deltal, double *deltar)
{
    int i;
    double *vnikx = malloc(4 * sizeof(double));
    double *work = malloc(16 * sizeof(double));
    int mflag = 0;

    int lenxk = n + 4;
    // # initialize output vectors
    for (i = 1; i < n + 1; i++)
    {
        y[i - 1] = 0.0;
        hs0[i - 1] = 0.0;
        hs1[i - 1] = 0.0;
        hs2[i - 1] = 0.0;
        hs3[i - 1] = 0.0;
    }

    int ileft = 1;
    double eps = 1e-9; // eps = .1d-9

    for (i = 1; i < k + 1; i++)
    {
        ileft = interv(xknot, n + 1, x[i - 1], 0, 0, ileft, &mflag);
        if (mflag == 1)
        {
            if (x[i - 1] <= (xknot[ileft - 1] + eps))
            {
                ileft -= 1;
            }
            else
            {
                free(vnikx);
                free(work);
                return;
            }
        }
        bsplvd(xknot, lenxk, 4, x[i - 1], ileft, work, vnikx, 1, j, deltal, deltar);

        *j = ileft - 4 + 1;
        y[*j - 1] = y[*j - 1] + pow(w[i - 1], 2) * z[i - 1] * vnikx[1 - 1];
        hs0[*j - 1] = hs0[*j - 1] + pow(w[i - 1], 2) * pow(vnikx[1 - 1], 2);
        hs1[*j - 1] = hs1[*j - 1] + pow(w[i - 1], 2) * vnikx[1 - 1] * vnikx[2 - 1];
        hs2[*j - 1] = hs2[*j - 1] + pow(w[i - 1], 2) * vnikx[1 - 1] * vnikx[3 - 1];
        hs3[*j - 1] = hs3[*j - 1] + pow(w[i - 1], 2) * vnikx[1 - 1] * vnikx[4 - 1];

        *j = ileft - 4 + 2;
        y[*j - 1] = y[*j - 1] + pow(w[i - 1], 2) * z[i - 1] * vnikx[2 - 1];
        hs0[*j - 1] = hs0[*j - 1] + pow(w[i - 1], 2) * pow(vnikx[2 - 1], 2);
        hs1[*j - 1] = hs1[*j - 1] + pow(w[i - 1], 2) * vnikx[2 - 1] * vnikx[3 - 1];
        hs2[*j - 1] = hs2[*j - 1] + pow(w[i - 1], 2) * vnikx[2 - 1] * vnikx[4 - 1];

        *j = ileft - 4 + 3;
        y[*j - 1] = y[*j - 1] + pow(w[i - 1], 2) * z[i - 1] * vnikx[3 - 1];
        hs0[*j - 1] = hs0[*j - 1] + pow(w[i - 1], 2) * pow(vnikx[3 - 1], 2);
        hs1[*j - 1] = hs1[*j - 1] + pow(w[i - 1], 2) * vnikx[3 - 1] * vnikx[4 - 1];

        *j = ileft - 4 + 4;
        y[*j - 1] = y[*j - 1] + pow(w[i - 1], 2) * z[i - 1] * vnikx[4 - 1];
        hs0[*j - 1] = hs0[*j - 1] + pow(w[i - 1], 2) * pow(vnikx[4 - 1], 2);
    }

    free(vnikx);
    free(work);
    return;
}

static double ddot(int n, double *dx, int incx, double *dy, int incy)
{
    // http://for2html.sourceforge.net/BLAS/ddot.html

    int i;
    double ddot = 0.0;
    double dtemp = 0.0;
    if (n <= 0)
        return ddot;

    if (incx == 1 && incy == 1)
    {
        int m, mp1;
        // 20
        m = n % 5;
        if (m != 0)
        {
            for (i = 1; i < m + 1; i++)
                dtemp = dtemp + dx[i - 1] * dy[i - 1];

            if (n < 5)
            {
                ddot = dtemp;
                return ddot;
            }
        }

        // 40
        mp1 = m + 1;
        for (i = mp1; i < n + 1; i += 5)
        {
            dtemp = (dtemp + dx[i - 1] * dy[i - 1] + dx[i + 1 - 1] * dy[i + 1 - 1] + dx[i + 2 - 1] * dy[i + 2 - 1] + dx[i + 3 - 1] * dy[i + 3 - 1] + dx[i + 4 - 1] * dy[i + 4 - 1]);
        }

        ddot = dtemp;
        return ddot;
    }

    int ix = 1; // = 1
    int iy = 1; // = 1
    if (incx < 0)
        ix = (-n + 1) * incx; // + 1
    if (incy < 0)
        iy = (-n + 1) * incy;
    for (i = 1; i < n + 1; i++)
    {
        dtemp = dtemp + dx[ix - 1] * dy[iy - 1];
        ix = ix + incx;
        iy = iy + incy;
    }

    ddot = dtemp;
    return ddot;
}

static void dpbfa(double *abd, int ld4, int n, int m, int info)
{
    // dpbfa factors a double precision symmetric positive definite
    int j, k, ik, jk, mu;
    double s;

    for (j = 1; j < n + 1; j++)
    {
        info = j;
        s = 0.0;
        ik = m + 1;
        jk = MAX(j - m, 1);
        mu = MAX(m + 2 - j, 1);
        if (!(m < mu))
            // go to 20
            for (k = mu; k < m + 1; k++)
            {
                double *x = abd + convert_2d_to_1d_idx(ld4, ik - 1, jk - 1);
                double *y = abd + convert_2d_to_1d_idx(ld4, mu - 1, j - 1);
                double t = abd[convert_2d_to_1d_idx(ld4, k - 1, j - 1)] - ddot(k - mu, x, 1, y, 1);
                t = t / abd[convert_2d_to_1d_idx(ld4, m + 1 - 1, jk - 1)];
                abd[convert_2d_to_1d_idx(ld4, k - 1, j - 1)] = t;
                s = s + t * t;
                ik = ik - 1;
                jk = jk + 1;
            }

        // 20 continue
        s = abd[convert_2d_to_1d_idx(ld4, m + 1 - 1, j - 1)] - s;
        if (s <= 0.0)
            break;

        abd[convert_2d_to_1d_idx(ld4, m + 1 - 1, j - 1)] = sqrt(s);
    }

    return;
}

static void daxpy(int n, double da, double *dx, int incx, double *dy, int incy)
{
    // http://for2html.sourceforge.net/BLAS/daxpy.html

    int i, m, mp1;

    if (n <= 0)
        return;
    if (da == 0.0)
        return;
    if (incx == 1 && incy == 1)
    {
        m = n % 4;
        if (m != 0)
        {
            for (i = 1; i < m + 1; i++)
                dy[i - 1] = dy[i - 1] + da * dx[i - 1];

            if (n < 4)
                return;
        }

        // 40
        mp1 = m + 1;
        for (i = mp1; i < n + 1; i += 4)
        {
            dy[i - 1] = dy[i - 1] + da * dx[i - 1];
            dy[i + 1 - 1] = dy[i + 1 - 1] + da * dx[i + 1 - 1];
            dy[i + 2 - 1] = dy[i + 2 - 1] + da * dx[i + 2 - 1];
            dy[i + 3 - 1] = dy[i + 3 - 1] + da * dx[i + 3 - 1];
        }
    }

    int ix = 1;
    int iy = 1;
    if (incx < 0)
        ix = (-n + 1) * incx + 1;
    if (incy < 0)
        iy = (-n + 1) * incy + 1;
    for (i = 1; i < n + 1; i++)
    {
        dy[iy - 1] = dy[iy - 1] + da * dx[ix - 1];
        ix = ix + incx;
        iy = iy + incy;
    }
}

static void dpbsl(double *abd, int lda, int n, int m, double *b)
{
    int k, kb, lm, la, lb;
    double *x;
    double t;

    for (k = 1; k < n + 1; k++)
    {
        lm = MIN(k - 1, m);
        la = m + 1 - lm;
        lb = k - lm;

        x = abd + convert_2d_to_1d_idx(lda, la - 1, k - 1);
        t = ddot(lm, x, 1, b + lb - 1, 1);
        b[k - 1] = (b[k - 1] - t) / abd[convert_2d_to_1d_idx(lda, m + 1 - 1, k - 1)];
    }

    // solve r*x = y
    for (kb = 1; kb < n + 1; kb++)
    {
        k = n + 1 - kb;
        lm = MIN(k - 1, m);
        la = m + 1 - lm;
        lb = k - lm;
        b[k - 1] = b[k - 1] / abd[convert_2d_to_1d_idx(lda, m + 1 - 1, k - 1)];
        t = -1.0 * b[k - 1];

        x = abd + convert_2d_to_1d_idx(lda, la - 1, k - 1);
        daxpy(lm, t, x, 1, b + lb - 1, 1);
    }
}

static double bvalue(double *t, double *bcoef, int n, int k, double x, int jderiv)
{
    int kmax = 20;
    double *aj = malloc(kmax * sizeof(double));
    double *dm = malloc(kmax * sizeof(double));
    double *dp = malloc(kmax * sizeof(double));
    int mflag, j;

    // initialize
    int i = 1;
    double bvalue = 0.0;

    if (jderiv >= k)
        return bvalue;

    if ((x != t[n + 1 - 1]) || (t[n + 1 - 1] != t[n + k - 1]))
    {
        i = interv(t, n + k, x, 0, 0, i, &mflag);
        if (mflag != 0)
            return bvalue;
    }
    else
    {
        i = n;
    }

    int km1 = k - 1;
    if (km1 <= 0)
    {
        bvalue = bcoef[i - 1];
        return bvalue;
    }

    int jcmin = 1;
    int imk = i - k;
    if (imk >= 0)
    {
        for (j = 1; j < km1 + 1; j++)
        {
            dm[j - 1] = x - t[i + 1 - j - 1];
        }
    }
    else
    {
        jcmin = 1 - imk;
        for (j = 1; j < i + 1; j++)
            dm[j - 1] = x - t[i + 1 - j - 1];

        for (j = i; j < km1 + 1; j++)
        {
            aj[k - j - 1] = 0.0;
            dm[j - 1] = dm[i - 1];
        }
    }

    int jcmax = k;
    int nmi = n - i;
    if (nmi >= 0)
    {
        for (j = 1; j < km1 + 1; j++)
            dp[j - 1] = t[i + j - 1] - x;
    }
    else
    {
        jcmax = k + nmi;
        for (j = 1; j < jcmax + 1; j++)
            dp[j - 1] = t[i + j - 1] - x;

        for (j = jcmax; j < km1 + 1; j++)
        {
            aj[j + 1 - 1] = 0.0;
            dp[j - 1] = dp[jcmax - 1];
        };
    }

    int jj, jc, ilo, kmj;
    for (jc = jcmin; jc < jcmax + 1; jc++)
        aj[jc - 1] = bcoef[imk + jc - 1];

    if (jderiv >= 1)
        for (j = 1; j < jderiv + 1; j++)
        {
            kmj = k - j;
            double fkmj = (double)kmj;
            ilo = kmj;
            for (jj = 1; jj < kmj + 1; jj++)
            {
                aj[jj - 1] = ((aj[jj + 1 - 1] - aj[jj - 1]) / (dm[ilo - 1] + dp[jj - 1])) * fkmj;
                ilo = ilo - 1;
            }
        }

    if (jderiv != km1)
    {
        int jdrvp1 = jderiv + 1;
        for (j = jdrvp1; j < km1 + 1; j++)
        {
            kmj = k - j;
            ilo = kmj;
            for (jj = 1; jj < kmj + 1; jj++)
            {
                aj[jj - 1] = (aj[jj + 1 - 1] * dm[ilo - 1] + aj[jj - 1] * dp[jj - 1]) / (dm[ilo - 1] + dp[jj - 1]);
                ilo = ilo - 1;
            }
        }
    }

    bvalue = aj[1 - 1];

    free(aj);
    free(dm);
    free(dp);
    return bvalue;
}

static void sinerp(
    double *abd,
    int ld4,
    int nk,
    double *p1ip,
    double *p2ip,
    int ldnk,
    int flag)
{
    double *wjm3 = malloc(3 * sizeof(double));
    double *wjm2 = malloc(2 * sizeof(double));
    double *wjm1 = malloc(1 * sizeof(double));

    double c0;
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0;

    // Pass 1
    wjm3[1 - 1] = 0.0;
    wjm3[2 - 1] = 0.0;
    wjm3[3 - 1] = 0.0;
    wjm2[1 - 1] = 0.0;
    wjm2[2 - 1] = 0.0;
    wjm1[1 - 1] = 0.0;

    int i, j, k;
    for (i = 1; i < nk + 1; i++)
    {
        j = nk - i + 1;
        c0 = 1.0 / abd[convert_2d_to_1d_idx(ld4, 4 - 1, j - 1)];
        if (j <= nk - 3)
        {
            c1 = abd[convert_2d_to_1d_idx(ld4, 1 - 1, j + 3 - 1)] * c0;
            c2 = abd[convert_2d_to_1d_idx(ld4, 2 - 1, j + 2 - 1)] * c0;
            c3 = abd[convert_2d_to_1d_idx(ld4, 3 - 1, j + 1 - 1)] * c0;
        }
        else if (j == nk - 2)
        {
            c1 = 0.0;
            c2 = abd[convert_2d_to_1d_idx(ld4, 2 - 1, j + 2 - 1)] * c0;
            c3 = abd[convert_2d_to_1d_idx(ld4, 3 - 1, j + 1 - 1)] * c0;
        }
        else if (j == nk - 1)
        {
            c1 = 0.0;
            c2 = 0.0;
            c3 = abd[convert_2d_to_1d_idx(ld4, 3 - 1, j + 1 - 1)] * c0;
        }
        else if (j == nk)
        {
            c1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;
        }

        p1ip[convert_2d_to_1d_idx(ld4, 1 - 1, j - 1)] = 0.0 - (c1 * wjm3[1 - 1] + c2 * wjm3[2 - 1] + c3 * wjm3[3 - 1]);
        p1ip[convert_2d_to_1d_idx(ld4, 2 - 1, j - 1)] = 0.0 - (c1 * wjm3[2 - 1] + c2 * wjm2[1 - 1] + c3 * wjm2[2 - 1]);
        p1ip[convert_2d_to_1d_idx(ld4, 3 - 1, j - 1)] = 0.0 - (c1 * wjm3[3 - 1] + c2 * wjm2[2 - 1] + c3 * wjm1[1 - 1]);
        p1ip[convert_2d_to_1d_idx(ld4, 4 - 1, j - 1)] = (pow(c0, 2) + pow(c1, 2) * wjm3[1 - 1] + 2.0 * c1 * c2 * wjm3[2 - 1] + 2.0 * c1 * c3 * wjm3[3 - 1] + pow(c2, 2) * wjm2[1 - 1] + 2.0 * c2 * c3 * wjm2[2 - 1] + pow(c3, 2) * wjm1[1 - 1]);
        wjm3[1 - 1] = wjm2[1 - 1];
        wjm3[2 - 1] = wjm2[2 - 1];
        wjm3[3 - 1] = p1ip[convert_2d_to_1d_idx(ld4, 2 - 1, j - 1)];
        wjm2[1 - 1] = wjm1[1 - 1];
        wjm2[2 - 1] = p1ip[convert_2d_to_1d_idx(ld4, 3 - 1, j - 1)];
        wjm1[1 - 1] = p1ip[convert_2d_to_1d_idx(ld4, 4 - 1, j - 1)];

        if (flag != 0)
        { // ____ Pass 2 _____  Compute p2ip  [never from R's code!]
            for (i = 1; i < nk + 1; i++)
                j = nk - i + 1;
            for (k = 1; k < 4 + 1; k++)
            {
                if (j + k - 1 > nk)
                    break;
                p2ip[convert_2d_to_1d_idx(ldnk, j - 1, j + k - 1 - 1)] = p1ip[convert_2d_to_1d_idx(ld4, 5 - k - 1, j - 1)];
            }

            for (i = 1; i < nk + 1; i++)
            {
                j = nk - i + 1;

                if (j - 4 >= 1)
                    for (k = j - 4; k > 1 + 1; k -= 1)
                    {
                        c0 = 1.0 / abd[convert_2d_to_1d_idx(ld4, 4 - 1, k - 1)];
                        c1 = abd[convert_2d_to_1d_idx(ld4, 1 - 1, k + 3 - 1)] * c0;
                        c2 = abd[convert_2d_to_1d_idx(ld4, 2 - 1, k + 2 - 1)] * c0;
                        c3 = abd[convert_2d_to_1d_idx(ld4, 3 - 1, k + 1 - 1)] * c0;
                        p2ip[convert_2d_to_1d_idx(ldnk, k - 1, j - 1)] = 0.0 - (c1 * p2ip[convert_2d_to_1d_idx(ldnk, k + 3 - 1, j - 1)] + c2 * p2ip[convert_2d_to_1d_idx(ldnk, k + 2 - 1, j - 1)] + c3 * p2ip[convert_2d_to_1d_idx(ldnk, k + 1 - 1, j - 1)]);
                    }
            }
        }
    }
    free(wjm3);
    free(wjm2);
    free(wjm1);
    return;
}

static void sslvrg(int *penalt, int *dofoff,
                   double *xs, double *ys, double *ws, double *ssw, int n,
                   double *knot, int nk,
                   double *coef, double *sz, double *lev, double *crit,
                   int *icrit, double *lspar, double *xwy,
                   double *hs0, double *hs1, double *hs2, double *hs3,
                   double *sg0, double *sg1, double *sg2, double *sg3, double *abd,
                   double *p1ip, double *p2ip, int *ld4, int *ldnk, int *ier,
                   int *j, double *deltal, double *deltar)
{
    // p1ip.shape = (ld4_ptr[0], nk)

    // local variables
    int mflag = 0;
    double *vnikx = malloc(4 * sizeof(double));
    double *work = malloc(16 * sizeof(double));

    int lenkno = nk + 4;
    int ileft = 1;
    double eps = 1e-11; // 1d-11

    // compute the coefficients coef() of estimated smooth
    int i;
    double xv;
    for (i = 1; i < nk + 1; i++)
    {
        coef[i - 1] = xwy[i - 1];
        abd[convert_2d_to_1d_idx(*ld4, 4 - 1, i - 1)] = (hs0[i - 1] + lspar[0] * sg0[i - 1]);
    }

    for (i = 1; i < nk - 1 + 1; i++)
    {
        abd[convert_2d_to_1d_idx(*ld4, 3 - 1, i + 1 - 1)] = (hs1[i - 1] + lspar[0] * sg1[i - 1]);
    }

    for (i = 1; i < nk - 2 + 1; i++)
    {
        abd[convert_2d_to_1d_idx(*ld4, 2 - 1, i + 2 - 1)] = (hs2[i - 1] + lspar[0] * sg2[i - 1]);
    }

    for (i = 1; i < nk - 3 + 1; i++)
    {
        abd[convert_2d_to_1d_idx(*ld4, 1 - 1, i + 3 - 1)] = (hs3[i - 1] + lspar[0] * sg3[i - 1]);
    }

    // factorize banded matrix abd (into upper triangular)
    dpbfa(abd, *ld4, nk, 3, *ier);
    if (*ier != 0) // matrix could not be factorized -> ier := info
        return;

    // solve linear system (from factorized abd):
    dpbsl(abd, *ld4, nk, 3, coef);

    // Value of smooth at the data points
    for (i = 1; i < n + 1; i++)
    {
        xv = xs[i - 1];
        sz[i - 1] = bvalue(knot, coef, nk, 4, xv, 0);
    }

    // Compute the criterion function if requested (icrit > 0)
    if (*icrit >= 1)
    {
        // Ordinary or Generalized CV or "df match"
        // Get Leverages First
        sinerp(abd, *ld4, nk, p1ip, p2ip, *ldnk, 0);
        for (i = 1; i < n + 1; i++)
        {
            xv = xs[i - 1];
            ileft = interv(knot, nk + 1, xv, 0, 0, ileft, &mflag);
            if (mflag == -1)
            {
                ileft = 4;
                xv = knot[4 - 1] + eps;
            }
            else if (mflag == 1)
            {
                ileft = nk;
                xv = knot[nk + 1 - 1] - eps;
            }

            int j_idx = ileft - 3;

            bsplvd(
                knot,
                lenkno,
                4,
                xv,
                ileft,
                work,
                vnikx,
                1,
                j,
                deltal,
                deltar);
            double b0 = vnikx[1 - 1];
            double b1 = vnikx[2 - 1];
            double b2 = vnikx[3 - 1];
            double b3 = vnikx[4 - 1];
            lev[i - 1] = (p1ip[convert_2d_to_1d_idx(*ld4, 4 - 1, j_idx - 1)] * pow(b0, 2) + 2.0 * p1ip[convert_2d_to_1d_idx(*ld4, 3 - 1, j_idx - 1)] * b0 * b1 + 2.0 * p1ip[convert_2d_to_1d_idx(*ld4, 2 - 1, j_idx - 1)] * b0 * b2 + 2.0 * p1ip[convert_2d_to_1d_idx(*ld4, 1 - 1, j_idx - 1)] * b0 * b3 + p1ip[convert_2d_to_1d_idx(*ld4, 4 - 1, j_idx + 1 - 1)] * pow(b1, 2) + 2.0 * p1ip[convert_2d_to_1d_idx(*ld4, 3 - 1, j_idx + 1 - 1)] * b1 * b2 + 2.0 * p1ip[convert_2d_to_1d_idx(*ld4, 2 - 1, j_idx + 1 - 1)] * b1 * b3 + p1ip[convert_2d_to_1d_idx(*ld4, 4 - 1, j_idx + 2 - 1)] * pow(b2, 2) + 2.0 * p1ip[convert_2d_to_1d_idx(*ld4, 3 - 1, j_idx + 2 - 1)] * b2 * b3 + p1ip[convert_2d_to_1d_idx(*ld4, 4 - 1, j_idx + 3 - 1)] * pow(b3, 2)) * pow(ws[i - 1], 2);
        }

        // Evaluate Criterion
        double df = 0.0;
        if (*icrit == 1) // ! Generalized CV --------------------
        {
            double rss = *ssw;
            double sumw = 0.0;

            for (i = 1; i < n + 1; i++)
            {
                rss = rss + pow(((ys[i - 1] - sz[i - 1]) * ws[i - 1]), 2);
                df = df + lev[i - 1];
                sumw = sumw + pow(ws[i - 1], 2);
            }

            *crit = (rss / sumw) / pow((1.0 - (*dofoff + *penalt * df) / sumw), 2);
        }
        else if (*icrit == 2) // ! Ordinary CV ------------------
        {
            *crit = 0.0;
            for (i = 1; i < n + 1; i++)
            {
                *crit = (*crit + pow((((ys[i - 1] - sz[i - 1]) * ws[i - 1]) / (1 - lev[i - 1])), 2));
            }
            *crit = *crit / n;
        }
        else
        { // df := sum( lev[i] )
            for (i = 1; i < n + 1; i++)
                df = df + lev[i - 1];

            if (*icrit == 3) // ! df matching --------------------
            {
                *crit = 3 + pow(*dofoff - df, 2);
            }
            else
            { // ! if(icrit .eq. 4) then df - dofoff (=> zero finding)
                *crit = df - *dofoff;
            }
        }
    }

    free(vnikx);
    free(work);
    return;
}

static void sbart(int *penalt, int *dofoff,
                  double *xs, double *ys, double *ws, double *ssw,
                  int n, double *knot, int nk, double *coef,
                  double *sz, double *lev, double *crit,
                  int *icrit, double *spar, int *ispar, int *iter,
                  double *lspar, double *uspar, double *tol, double *eps, double *Ratio,
                  int *isetup,
                  double *xwy, double *hs0, double *hs1, double *hs2,
                  double *hs3, double *sg0, double *sg1, double *sg2,
                  double *sg3, double *abd, double *p1ip, double *p2ip,
                  int *ld4, int *ldnk, int *ier)
{
    // "Correct" ./sslvrg.f (line 129):   crit = 3 + (dofoff-df)**2
#define CRIT(FX) (*icrit == 3 ? FX - 3. : FX)
    /* cancellation in (3 + eps) - 3, but still...informative */

#define BIG_f (1e100)

    /* c_Gold is the squared inverse of the golden ratio */
    static const double c_Gold = 0.381966011250105151795413165634;
    /* == (3. - sqrt(5.)) / 2. */

    /* Local variables */
    static double ratio; /* must be static (not needed in R) */

    double a, b, d, e, p, q, r, u, v, w, x;
    double ax, fu, fv, fw, fx, bx, xm;
    double tol1, tol2;

    int i, maxit;
    // uint8_t Fparabol = FALSE, tracing = (*ispar < 0),
    uint8_t spar_is_lambda = 0;

    /* unnecessary initializations to keep  -Wall happy */
    d = 0.;
    fu = 0.;
    u = 0.;
    // never computed if(spar_is_lambda)
    ratio = 1.;

    int j = 1;
    int jmax = 20;
    double *deltal = malloc(jmax * sizeof(double));
    double *deltar = malloc(jmax * sizeof(double));

    /*  Compute SIGMA, X' W X, X' W z, trace ratio, s0, s1.
        SIGMA	-> sg0,sg1,sg2,sg3   -- via sgram() in ./sgram.f
        X' W X	-> hs0,hs1,hs2,hs3   \
        X' W Z	-> xwy               _\ via stxwx() in ./stxwx.f
    */

    /* trevor fixed this 4/19/88
     * Note: sbart, i.e. stxwx() and sslvrg() {mostly, not always!}, use
     *	 the square of the weights; the following rectifies that */
    for (i = 0; i < n; ++i)
        if (ws[i] > 0.)
            ws[i] = sqrt(ws[i]);

    if (*isetup < 0)
        spar_is_lambda = 1;
    else if (*isetup != 1)
    { // 0 or 2
        /* SIGMA[i,j] := Int  B''(i,t) B''(j,t) dt  {B(k,.) = k-th B-spline} */

        sgram(sg0, sg1, sg2, sg3, knot, nk, &j, deltal, deltar);
        stxwx(xs, ys, ws, n,
              knot, nk,
              xwy,
              hs0, hs1, hs2, hs3,
              &j, deltal, deltar);
        spar_is_lambda = (*isetup == 2);
        if (!spar_is_lambda)
        {
            /* Compute ratio :=  tr(X' W X) / tr(SIGMA) */
            double t1 = 0., t2 = 0.;
            for (i = 3 - 1; i < (nk - 3); ++i)
            {
                t1 += hs0[i];
                t2 += sg0[i];
            }
            ratio = t1 / t2;
        }
        *isetup = 1;
    }
/*     Compute estimate */

// Compute SSPLINE(SPAR), assign result to *crit (and the auxil.variables)
#define SSPLINE_COMP(_SPAR_)                                         \
    *lspar = spar_is_lambda ? _SPAR_                                 \
                            : ratio * pow(16.0, (_SPAR_)*6.0 - 2.0); \
    sslvrg(penalt, dofoff, xs, ys, ws, ssw, n,                       \
           knot, nk,                                                 \
           coef, sz, lev, crit, icrit, lspar, xwy,                   \
           hs0, hs1, hs2, hs3,                                       \
           sg0, sg1, sg2, sg3, abd,                                  \
           p1ip, p2ip, ld4, ldnk, ier, &j, deltal, deltar)

    if (*ispar == 1)
    { /* Value of spar supplied */
        SSPLINE_COMP(*spar);
        /* got through check 2 */
        *Ratio = ratio;

        free(deltal);
        free(deltar);
        return;
    }

    /* ELSE ---- spar not supplied --> compute it ! ---------------------------
     */
    ax = *lspar;
    bx = *uspar;

    /* initialization */

    maxit = *iter;
    *iter = 0;
    a = ax;
    b = bx;
    v = a + c_Gold * (b - a);
    w = v;
    x = v;
    e = 0.;
    SSPLINE_COMP(x);
    fx = *crit;
    fv = fx;
    fw = fx;

    /* main loop
       --------- */
    while (*ier == 0)
    { /* L20: */
        xm = (a + b) * .5;
        tol1 = *eps * fabs(x) + *tol / 3.;
        tol2 = tol1 * 2.;
        ++(*iter);

        /* Check the (somewhat peculiar) stopping criterion: note that
           the RHS is negative as long as the interval [a,b] is not small:*/
        if (fabs(x - xm) <= tol2 - (b - a) * .5 || *iter > maxit)
            goto L_End;

        /* is golden-section necessary */

        if (fabs(e) <= tol1 ||
            /*  if had Inf then go to golden-section */
            fx >= BIG_f || fv >= BIG_f || fw >= BIG_f)
            goto L_GoldenSect;

        r = (x - w) * (fx - fv);
        q = (x - v) * (fx - fw);
        p = (x - v) * q - (x - w) * r;
        q = (q - r) * 2.;
        if (q > 0.)
            p = -p;
        q = fabs(q);
        r = e;
        e = d;

        /* is parabola acceptable?  Otherwise do golden-section */

        if (fabs(p) >= fabs(.5 * q * r) ||
            q == 0.)
            /* above line added by BDR;
             * [the abs(.) >= abs() = 0 should have branched..]
             * in FTN: COMMON above ensures q is NOT a register variable */

            goto L_GoldenSect;

        if (p <= q * (a - x) ||
            p >= q * (b - x))
            goto L_GoldenSect;

        /* Parabolic Interpolation step */
        d = p / q;
        u = x + d;

        /* f must not be evaluated too close to ax or bx */
        if (u - a < tol2 ||
            b - u < tol2)
            d = fsign(tol1, xm - x);

        goto L50;
        /*------*/

    L_GoldenSect: /* a golden-section step */
        if (x >= xm)
            e = a - x;
        else /* x < xm*/
            e = b - x;
        d = c_Gold * e;

    L50:
        u = x + ((fabs(d) >= tol1) ? d : fsign(tol1, d));
        /*  tol1 check : f must not be evaluated too close to x */

        SSPLINE_COMP(u);
        fu = *crit;
        // if (!R_FINITE(fu))
        // {
        //     fu = 2. * BIG_f;
        // }

        /*  update  a, b, v, w, and x */

        if (fu <= fx)
        {
            if (u >= x)
                a = x;
            else
                b = x;

            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            if (u < x)
                a = u;
            else
                b = u;

            if (fu <= fw || w == x)
            { /* L70: */
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            { /* L80: */
                v = u;
                fv = fu;
            }
        }
    } /* end main loop -- goto L20; */

L_End:
    *Ratio = ratio;
    *spar = x;
    *crit = fx;

    free(deltal);
    free(deltar);
    return;
}

static void rbart(
    int *penalt,
    int *dofoff,
    double *xs,
    double *ys,
    double *ws,
    double *ssw,
    int n,
    double *knot,
    int nk,
    double *coef,
    double *sz,
    double *lev,
    double *crit,
    int *iparms,
    double *spar,
    double *lspar,
    double *uspar,
    double *tol,
    double *eps,
    double *Ratio,
    double *scrtch,
    int ld4,
    int ldnk,
    int *ier)
{
    int isetup = (iparms[3] == 1) ? 2 : 0;
    sbart(
        penalt,
        dofoff,
        xs,
        ys,
        ws,
        ssw,
        n,
        knot,
        nk,
        coef,
        sz,
        lev,
        crit,
        &iparms[0], // icrit
        spar,
        &iparms[1], // ispar
        &iparms[2], // _iter
        lspar,
        uspar,
        tol,
        eps,
        Ratio,
        &isetup,
        scrtch,                           // xwy
        scrtch + nk,                      // hs0
        scrtch + (2 * nk),                // hs1
        scrtch + (3 * nk),                // hs2
        scrtch + (4 * nk),                // hs3
        scrtch + (5 * nk),                // sg0
        scrtch + (6 * nk),                // sg1
        scrtch + (7 * nk),                // sg2
        scrtch + (8 * nk),                // sg3
        scrtch + (9 * nk),                // abd
        scrtch + (9 * nk + ld4 * nk),     // p1ip
        scrtch + (9 * nk + 2 * ld4 * nk), // p2ip
        &ld4,
        &ldnk,
        ier);
    return;
}

static double *bvalus(int n, double *knot, double *coef, int nk, double *x, double *s, int order)
{
    int i;
    for (i = 0; i < n; i++)
        s[i] = bvalue(knot, coef, nk, 4, x[i], order);

    return s;
}

void free_smooth_spline_model(smooth_spline_model_t *model)
{
    // free(model->x);          // Shouldn't free sm_spline_coef inputs
    free(model->y);
    // free(model->w);          // Shouldn't free sm_spline_coef inputs
    free(model->yin);
    // free(model->data_x);     // Shouldn't free sm_spline_coef inputs
    // free(model->data_y);     // Shouldn't free sm_spline_coef inputs
    // free(model->data_w);     // Shouldn't free sm_spline_coef inputs
    free(model->fit_knot);

    // free fit_t pointers
    fit_t *fit = model->fit;
    free(fit->coef);
    free(fit->lev);
    free(fit->parms);
    free(fit->crit);
    free(fit->iparms);
    free(fit->ier);
    free(fit->scrtch);

    free(fit);
    free(model);
}

smooth_spline_model_t *sm_spline_coef(int n, double *x, double *y, int w_len, double *w, double spar)
{
    int i;
    double tol = 1e-6 * iqr(n, x);

    contr_sp_t *contr_sp = malloc(sizeof(contr_sp_t));
    contr_sp->low = -1.5;
    contr_sp->high = 1.5;
    contr_sp->tol = 1e-4;
    contr_sp->eps = 2e-8;
    contr_sp->maxit = 500;
    contr_sp->trace = 0;

    uint8_t no_wgts = (w_len == 0);
    if (no_wgts)
    {
        w = malloc(n * sizeof(double));
        for (i = 0; i < n; i++)
        {
            w[i] = 1;
        }
    }
    else
    {
        double w_positive_sum = 0;
        double w_sum = 0;
        for (i = 0; i < n; i++)
        {
            if (w[i] > 0)
            {
                w_positive_sum += w[i];
            }
            w_sum += w[i];
        }
        for (i = 0; i < n; i++)
        {
            w[i] = (w[i] * w_positive_sum) / w_sum;
        }
    }

    // Assume all x unique
    double *ux = x;
    int nx = n;
    uint8_t ox = 1;
    // double tmp[n][3];
    // for (i = 0; i < n; i++)
    // {
    //     tmp[i][0] = w[i];
    //     tmp[i][1] = w[i] * y[i];
    //     tmp[i][2] = w[i] * (y[i] * y[i]);
    // }
    double *wbar = w;

    double *ybar = malloc(n * sizeof(double));
    uint8_t wbar_all_positive = 1;
    for (i = 0; i < n; i++)
    {
        if (wbar[i] <= 0)
        {
            wbar_all_positive = 0;
            break;
        }
    }

    if (wbar_all_positive)
    {
        for (i = 0; i < n; i++)
        {
            ybar[i] = (w[i] * y[i]) / wbar[i]; // bro this is just y[i] tbh
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            ybar[i] = w[i] * y[i];
        }
    }

    double yssw = 0;
    for (i = 0; i < n; i++)
    {
        yssw += (w[i] * (y[i] * y[i])) - (wbar[i] * (ybar[i] * ybar[i]));
    }

    double r_ux = ux[nx - 1] - ux[0];
    double *xbar = malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        xbar[i] = (ux[i] - ux[0]) / r_ux;
    }

    // all_knots = False by default
    int nknots = _nknots_smspl(nx);
    double *mid = malloc(nknots * sizeof(double));
    int step = (nx) / nknots;
    int xbar_idx = 0;
    for (i = 0; i < nknots; i++)
    {
        mid[i] = xbar[xbar_idx];
        xbar_idx += step;
    }
    int knot_size = 3 + nknots + 3;
    double *knot = malloc(knot_size * sizeof(double));
    for (i = 0; i < 3; i++)
    {
        knot[i] = xbar[0];
    }
    for (i = 3; i < nknots + 3; i++)
    {
        knot[i] = mid[i - 3];
    }
    free(mid);
    for (i = nknots + 3; i < nknots + 6; i++)
    {
        knot[i] = xbar[nx - 1];
    }

    int nk = nknots + 2;

    // spar stuff
    int spar_is_lambda = 0; // 0
    int ispar = 1;

    // icrit stuff
    int icrit = 1;
    int df_offset = 0;
    int dofoff = df_offset;
    double *coef = malloc(nk * sizeof(double));
    double *ty = malloc(nx * sizeof(double));
    double *lev = malloc(nx * sizeof(double));
    // double *lspar = &contr_sp.low;
    // double *uspar = &contr_sp.high;
    // double *tol = &contr_sp.tol;
    // double *eps = &contr_sp.eps;
    double Ratio = -1.0;

    double *crit = calloc(1, sizeof(double));
    int *ier = calloc(1, sizeof(int));
    int *iparms = malloc(4 * sizeof(int));
    iparms[0] = icrit;
    iparms[1] = ispar;
    iparms[2] = contr_sp->maxit;
    iparms[3] = spar_is_lambda;
    double *scrtch = malloc((18 * nk + 1) * sizeof(double));

    int penalty = 1;
    double *ws = malloc(w_len * sizeof(double));
    for (i = 0; i < w_len; i++)
    {
        ws[i] = wbar[i];
    }

    rbart(
        &penalty,        // penalt
        &dofoff,         // dofoff
        xbar,            // xs
        ybar,            // ys
        ws,              // ws
        &yssw,           // ssw
        nx,              // n
        knot,            // knot
        nk,              // nk
        coef,            // coef
        ty,              // sz
        lev,             // lev
        crit,            // crit
        iparms,          // iparms
        &spar,           // spar
        &contr_sp->low,  // lspar
        &contr_sp->high, // uspar
        &contr_sp->tol,  // tol
        &contr_sp->eps,  // eps
        &Ratio,          // Ratio
        scrtch,          // scrtch
        4,               // ld4
        1,               // ldnk
        ier              // ier
    );
    free(xbar);
    free(ws);

    double *parms = malloc(5 * sizeof(double));
    parms[0] = contr_sp->low;
    parms[1] = contr_sp->high;
    parms[2] = contr_sp->tol;
    parms[3] = contr_sp->eps;
    parms[4] = Ratio;

    fit_t *fit = malloc(sizeof(fit_t));
    fit->coef = coef;     // double *coef;
    fit->ty = ty;         // double *ty;
    fit->lev = lev;       // double *lev;
    fit->spar = spar;     // double *spar;
    fit->parms = parms;   // double *parms;
    fit->crit = crit;     // double *crit;
    fit->iparms = iparms; // int *iparms;
    fit->ier = ier;       // int *ier;
    fit->scrtch = scrtch; // double *scrtch;

    double df = 0;
    for (i = 0; i < nx; i++)
    {
        df += fit->lev[i];
    }

    smooth_spline_model_t *model = malloc(sizeof(smooth_spline_model_t));
    if (fit->ier[0] > 0)
    {
        uint8_t sml = fit->spar < 0.5;
        if (spar_is_lambda || sml)
        {
            return model;
        }
        else
        {
            fit->ty = malloc(nx * sizeof(double));
            double y_mean = mean(nx, y);
            for (i = 0; i < n; i++)
            {
                fit->ty[i] = y_mean;
            }
            df = 1;
        }
    }

    // cv is false, not none
    // double *r = malloc(n * sizeof(double));
    // for (i = 0; i < n; i++)
    // {
    //     r[i] = y[i] - fit->ty[i];
    // }
    double *r_squared = malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        r_squared[i] = y[i] - fit->ty[i];
        r_squared[i] = r_squared[i] * r_squared[i];
    }
    double cv_crit_div = 1 - (df_offset + penalty * df) / (double)n;
    double cv_crit = weighted_avg(n, r_squared, w) / (cv_crit_div * cv_crit_div);
    free(r_squared);

    model->x = ux;               // double *x;          n
    model->y = fit->ty;          // double *y;          n
    model->w = wbar;             // double *w;          n
    model->yin = ybar;           // double *yin;        n
    model->tol = tol;            // double tol;
    model->data_x = x;           // double *data_x;     n
    model->data_y = y;           // double *data_y;     n
    model->data_w = w;           // double *data_w;     n
    model->no_weights = no_wgts; // uint8_t no_weights;
    model->lev = fit->lev;       // double *lev;        nx
    model->cv_crit = cv_crit;    // double cv_crit;
    model->df = df;              // double df;
    // fit
    model->fit = fit;
    model->fit_coef = fit->coef;
    model->crit = fit->crit;       // double* crit;
    model->spar = fit->spar;       // double spar;
    model->ratio = fit->parms[4];  // double ratio;
    model->lambda = fit->parms[0]; // double lambda;
    model->fit_knot = knot;
    model->fit_nk = nk;
    model->fit_min = ux[0];
    model->fit_range = r_ux;

    free(contr_sp);
    return model;
}

double *predict_smooth_spline(smooth_spline_model_t *model, double *x, int x_len, int deriv)
{
    int i;
    double *xs = malloc(x_len * sizeof(double));
    double *y = xs;
    int *extrap_left = malloc(x_len * sizeof(int));
    int *extrap_right = malloc(x_len * sizeof(int));
    int *extrap = malloc(x_len * sizeof(int));
    int *interp = malloc(x_len * sizeof(int));
    int n = 0;
    int any_extrap_left = 0;
    int any_extrap_right = 0;

    for (i = 0; i < x_len; i++)
    {
        xs[i] = (x[i] - model->fit_min) / model->fit_range; // x scaled to [0,1]
        if (xs[i] < 0)                                      // extrap_left[i] = (xs[i] < 0) ? 1 : 0;
        {
            extrap_left[i] = 1;
            any_extrap_left = 1;
        }
        if (xs[i] > 1)
        {
            extrap_right[i] = 1;
            any_extrap_right = 1;
        }

        extrap[i] = extrap_left[i] || extrap_right[i];
        interp[i] = extrap[i] ? 0 : 1;
        n += interp[i];
    }

    if (n > 0)
    {
        double *bvalus_xs = malloc(n * sizeof(double));
        int j = 0;
        for (i = 0; i < x_len; i++)
        {
            if (interp[i])
            {
                bvalus_xs[j] = xs[i];
                j += 1;
            }
        }
        double *bvalus_result = bvalus(
            (int)n,
            model->fit_knot,
            model->fit_coef,
            (int)model->fit_nk,
            bvalus_xs,
            malloc(n * sizeof(double)),
            (int)deriv);

        free(bvalus_xs);

        j = 0;
        for (i = 0; i < x_len; i++)
        {
            if (interp[i])
            {
                y[i] = bvalus_result[j];
                j += 1;
            }
        }

        free(interp);
        free(bvalus_result);
    }

    if (n < x_len)
    {
        double xrange[2] = {model->fit_min, model->fit_min + model->fit_range};
        double *end_model_fit = malloc(2 * sizeof(double));
        double *end_slopes = malloc(2 * sizeof(double));

        if (deriv == 0)
        {
            end_model_fit = predict_smooth_spline(model, xrange, 2, 0);
            end_slopes = predict_smooth_spline(model, xrange, 2, 1);
            end_slopes[0] *= model->fit_range;
            end_slopes[1] *= model->fit_range;

            if (any_extrap_left)
            {
                for (i = 0; i < x_len; i++)
                {
                    if (extrap_left[i])
                        y[i] = end_model_fit[0] + end_slopes[0] * xs[i];
                }
            }

            if (any_extrap_right)
            {
                for (i = 0; i < x_len; i++)
                {
                    if (extrap_right[i])
                        y[i] = end_model_fit[1] + end_slopes[1] * (xs[i] - 1);
                }
            }
        }
        else if (deriv == 1)
        {
            end_slopes = predict_smooth_spline(model, xrange, 2, 1);
            end_slopes[0] *= model->fit_range;
            end_slopes[1] *= model->fit_range;

            for (i = 0; i < x_len; i++)
            {
                if (extrap_left[i])
                    y[i] = end_slopes[0];
                if (extrap_right[i])
                    y[i] = end_slopes[1];
            }
        }
        else
        {
            for (i = 0; i < x_len; i++)
            {
                if (extrap[i])
                    y[i] = 0;
            }
        }
    }
    if (deriv > 0)
        for (i = 0; i < x_len; i++)
            y[i] /= pow(model->fit_range, deriv);

    free(extrap);
    free(extrap_left);
    free(extrap_right);
    return y;
}
