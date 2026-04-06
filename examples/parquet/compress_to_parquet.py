import polars as pl
import time

print("Loading 80,000 frames from CSV...")
start = time.time()
df = pl.read_csv("synthetic_traj.csv")
print(f"Loaded in {time.time() - start:.3f} seconds!")

print("Compressing to Apache Parquet...")
df.write_parquet("synthetic_traj.parquet")

import os
csv_size = os.path.getsize("synthetic_traj.csv") / (1024 * 1024)
pq_size = os.path.getsize("synthetic_traj.parquet") / (1024 * 1024)

print(f"CSV Size: {csv_size:.2f} MB")
print(f"Parquet Size: {pq_size:.2f} MB")
