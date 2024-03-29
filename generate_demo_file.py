import time
import random

st = time.time()
sr = 100
total_time = 600 # seconds

total_n = total_time * sr

ts = st

filepath = "data/flat.csv"
with open(filepath, "w") as f:
    f.write("timestamp,x,y,z\n")
    for i in range(total_n):
        x = 0
        y = 1
        z = 0
        row = f"{ts:.3f},{x:.6f},{y:.6f},{z:.6f}\n"
        ts += 1.0 / sr
        f.write(row)

filepath = "data/random.csv"
with open(filepath, "w") as f:
    f.write("timestamp,x,y,z\n")
    for i in range(total_n):
        x = random.random()
        y = random.random()
        z = random.random()
        row = f"{ts:.3f},{x:.6f},{y:.6f},{z:.6f}\n"
        ts += 1.0 / sr
        f.write(row)
