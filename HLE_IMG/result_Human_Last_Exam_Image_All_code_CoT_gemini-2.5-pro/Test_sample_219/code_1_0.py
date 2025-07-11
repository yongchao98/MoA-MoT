import numpy as np

# 1. Identify the timestamps of the synchronized bursts from the graph.
# The timestamps are estimated from the center of each vertical band of high activity.
burst_timestamps = [2.5, 9.0, 17.0, 25.0, 36.0, 48.0]
print(f"Identified burst timestamps (s): {burst_timestamps}")

# 2. Calculate the inter-burst intervals (IBIs).
inter_burst_intervals = []
for i in range(len(burst_timestamps) - 1):
    interval = burst_timestamps[i+1] - burst_timestamps[i]
    inter_burst_intervals.append(interval)

print(f"Calculated inter-burst intervals (s): {inter_burst_intervals}")

# 3. Calculate the average inter-burst interval.
average_ibi = np.mean(inter_burst_intervals)

# 4. Print the calculation and the result rounded to two significant figures.
calculation_str = " + ".join(map(str, inter_burst_intervals))
print(f"Average IBI calculation: ({calculation_str}) / {len(inter_burst_intervals)}")
print(f"Average inter-burst interval: {average_ibi:.2g} s")