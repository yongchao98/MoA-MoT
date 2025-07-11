import numpy as np

# Step 1: Identify the timestamps of the synchronized bursts from the graph.
# These are the approximate center times of the dense vertical firing bands.
burst_timestamps = [2, 12.5, 21.5, 35.5, 48.5]
print(f"Identified burst timestamps (s): {burst_timestamps}")
print("-" * 30)

# Step 2: Calculate the inter-burst intervals (IBIs).
# The IBI is the time difference between consecutive bursts.
intervals = np.diff(burst_timestamps)
print("Calculated Inter-Burst Intervals (s):")
for i, interval in enumerate(intervals):
    print(f"Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {interval} s")
print("-" * 30)

# Step 3: Calculate the average of the IBIs.
average_interval = np.mean(intervals)

# Step 4: Display the full calculation for the average.
interval_sum_str = " + ".join(map(str, intervals))
print("Average Inter-Burst Interval Calculation:")
print(f"({interval_sum_str}) / {len(intervals)} = {average_interval} s")
print("-" * 30)

# Step 5: Round the result to two significant figures and print.
# Using format specifier 'g' for significant figures.
average_interval_rounded = float(f"{average_interval:.2g}")
print(f"The average inter-burst interval is {average_interval} s.")
print(f"Rounded to two significant figures, the average inter-burst interval is {average_interval_rounded} s.")