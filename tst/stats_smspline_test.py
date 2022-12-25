import ctypes
import numpy as np
import pathlib
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from scipy.interpolate import interp1d
import sys

sys.path.insert(0, "../../..")
from pymims.helper import timestamp_sequence
from pymims.smooth_spline.smooth_spline import (
    predict_smooth_spline,
    smooth_spline_model,
)

# get inputs
spar = 0.6
over_t = np.array(
    [
        1485471758.4100000858306884765625000000000000000000000000000000000000000000000000000000000000000000000000000000,
        1485471758.4200000762939453125000000000000000000000000000000000000000000000000000000000000000000000000000000000,
        1485471758.4300000667572021484375000000000000000000000000000000000000000000000000000000000000000000000000000000,
        1485471758.4400000572204589843750000000000000000000000000000000000000000000000000000000000000000000000000000000,
        1485471758.4500000476837158203125000000000000000000000000000000000000000000000000000000000000000000000000000000,
    ]
)
over_values = np.array(
    [
        -2.7232404528720377534511953854234889149665832519531250000000000000000000000000000000000000000000000000,
        -2.9982694914466998525881535897497087717056274414062500000000000000000000000000000000000000000000000000,
        -3.2761937695231200251555492286570370197296142578125000000000000000000000000000000000000000000000000000,
        -3.5578875828916785017952406633412465453147888183593750000000000000000000000000000000000000000000000000,
        -3.8418727509985299839456729387165978550910949707031250000000000000000000000000000000000000000000000000,
    ]
)
weights = np.array(
    [
        0.8498633849249214167187460589047987014055252075195312500000000000000000000000000000000000000000000000,
        0.8498633849249214167187460589047987014055252075195312500000000000000000000000000000000000000000000000,
        0.8498633849249214167187460589047987014055252075195312500000000000000000000000000000000000000000000000,
        0.8498633849249214167187460589047987014055252075195312500000000000000000000000000000000000000000000000,
        1.6005464603003138890358059143181890249252319335937500000000000000000000000000000000000000000000000000,
    ]
)

# run R version
# r_smooth_spline = robjects.r["smooth.spline"]
# r_fitted = r_smooth_spline(
#     x=robjects.FloatVector(over_t),
#     y=robjects.FloatVector(over_values),
#     w=robjects.FloatVector(weights),
#     s=spar,
# )
# r_fitted = dict(zip(r_fitted.names, list(r_fitted)))

# run python version
# py_fitted = smooth_spline_model(x=over_t, y=over_values, w=weights, spar=spar)
py_fitted = smooth_spline_model(x=over_t, y=over_values, w=weights, spar=spar)["fit"]
py_results = predict_smooth_spline(py_fitted, over_t)["y"]

# Load the shared library into ctypes
# stats_smspline_lib = ctypes.CDLL(
#     "/Users/arytonhoi/Kode/mhealth/pymims/planning/C/spline/stats_smspline.so"
# )

# double_array_1 = ctypes.c_double * len(over_t)
# c_int = ctypes.c_int

# x_p = double_array_1(*over_t)
# y_p = double_array_1(*over_values)
# w_p = double_array_1(*weights)
# stats_smspline_lib.SplineCoef(
#     len(over_t), x_p, y_p, len(weights), w_p, ctypes.c_double(spar)
# )

# # compare
# expected = pd.read_csv(
#     "/Users/arytonhoi/Kode/mhealth/pymims/planning/C/spline/py_spline_eval.csv"
# )
# actual = pd.read_csv(
#     "/Users/arytonhoi/Kode/mhealth/pymims/planning/C/spline/spline_eval.csv"
# )

# # Mismatched elements: 293394 / 359996 (81.5%)
# # Max absolute difference: 6.66133815e-16
# # Max relative difference: 7.40517225e-11
# np.testing.assert_array_equal(
#     expected[expected.columns[0]].values, actual[actual.columns[1]].values[1:]
# )
