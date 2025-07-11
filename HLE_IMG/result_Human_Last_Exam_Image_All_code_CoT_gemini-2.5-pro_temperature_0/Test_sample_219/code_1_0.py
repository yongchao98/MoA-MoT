import numpy as np

# Step 1 & 2: Identify the timestamps of the synchronized bursts from the graph.
# The bursts are the dense vertical bands of activity.
# I estimate the center of each burst from the x-axis (Time (s)).
burst_timestamps = [2.5, 12.5, 24.5, 36.5, 48.5]
print(f"Identified burst timestamps (s): {burst_timestamps}")

# Step 3: Calculate the inter-burst intervals (IBIs).
# IBI is the time difference between consecutive bursts.
inter_burst_intervals = []
for i in range(len(burst_timestamps) - 1):
    interval = burst_timestamps[i+1] - burst_timestamps[i]
    inter_burst_intervals.append(interval)

print(f"Calculated inter-burst intervals (s): {inter_burst_intervals}")

# Step 4: Calculate the average IBI.
average_ibi = np.mean(inter_burst_intervals)

# Step 5: Format the output and print the final calculation.
# The final equation shows the sum of intervals divided by the number of intervals.
# The result is then formatted to two significant figures.
print("\nCalculating the average inter-burst interval:")
calculation_str = " + ".join(map(str, inter_burst_intervals))
print(f"Average = ({calculation_str}) / {len(inter_burst_intervals)}")
print(f"Average = {sum(inter_burst_intervals)} / {len(inter_burst_intervals)}")
print(f"Average = {average_ibi:.3f} s")

# Rounding to two significant figures
final_answer = float(f"{average_ibi:.1g}")
print(f"\nThe average inter-burst interval rounded to two significant figures is {final_answer} s.")