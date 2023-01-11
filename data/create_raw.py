import numpy as np
import pandas as pd

folder = "/Users/arytonhoi/Kode/mhealth/cmims/data/mims_unit/aditya_sleep"

timestamps = np.genfromtxt(f"{folder}/timestamps.csv", delimiter=",", dtype=np.int64)
timestamps = np.divide(timestamps, 10 ** 9)
x = np.genfromtxt(f"{folder}/x.csv", delimiter=",")
y = np.genfromtxt(f"{folder}/y.csv", delimiter=",")
z = np.genfromtxt(f"{folder}/z.csv", delimiter=",")

df = pd.DataFrame(
    {"timestamps": timestamps, "x": x, "y": y, "z": z},
    columns=["timestamps", "x", "y", "z"],
)

df.round(4)
df.to_csv(f"{folder}/raw.csv", index=False)
