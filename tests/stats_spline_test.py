import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import sys

sys.path.insert(0, "../../..")
from pymims.helper import timestamp_sequence
from pymims.spline import spline
from test.util import get_test_dataframe_pairs

mims_unit = importr("MIMSunit")

# Setup input dataframes
df_pairs = get_test_dataframe_pairs()
df_1 = df_pairs[0]["python"].copy()
df_2 = df_pairs[1]["python"].copy()

import ctypes
import pathlib

# get inputs
timestamps = df_1.index
float_timestamps = np.around(timestamps.values.astype(float) / (10 ** 9), 3)
oversampled_float_timestamps = timestamp_sequence(
    float_timestamps[0], float_timestamps[-1], step=1 / 100
)
values = df_1[df_1.columns[0]].values
x = float_timestamps
y = values
xout = oversampled_float_timestamps

# run R version
r_x = robjects.FloatVector(x)
r_y = robjects.FloatVector(y)
r_xout = robjects.FloatVector(xout)
r_spline = robjects.r["spline"]
# py_yout = np.array(r_spline(x=r_x, y=r_y, xout=r_xout, method="natural"))[1]
py_yout = np.array(r_spline(x=r_x, y=r_y, xout=r_xout, method="fmm"))[1]
np.savetxt("py_spline_eval.csv", py_yout, delimiter=",", fmt="%.15e")
# for o in py_yout:
#     print(f"{o:0.15f}")

# Load the shared library into ctypes
stats_spline_lib = ctypes.CDLL(
    "/Users/arytonhoi/Kode/mhealth/pymims/planning/C/spline/stats_spline.so"
)
# c_double_p = ctypes.POINTER(ctypes.c_double * 20)
# stats_spline_lib.spline_coef.restype = c_double_p

double_array_1 = ctypes.c_double * len(x)
double_array_2 = ctypes.c_double * len(xout)
x_p = double_array_1(*x)
y_p = double_array_1(*y)
xout_p = double_array_2(*xout)
# stats_spline_lib.Spline(len(x), x_p, len(y), y_p, len(xout), xout_p, 2)
stats_spline_lib.Spline(len(x), x_p, len(y), y_p, len(xout), xout_p, 3)

# compare
expected = pd.read_csv(
    "/Users/arytonhoi/Kode/mhealth/pymims/planning/C/spline/py_spline_eval.csv"
)
actual = pd.read_csv(
    "/Users/arytonhoi/Kode/mhealth/pymims/planning/C/spline/spline_eval.csv"
)

# Mismatched elements: 293394 / 359996 (81.5%)
# Max absolute difference: 6.66133815e-16
# Max relative difference: 7.40517225e-11
np.testing.assert_array_equal(
    expected[expected.columns[0]].values, actual[actual.columns[1]].values[1:]
)
