import numpy as np

# Step 1: Identify the timestamps of synchronized bursts from the graph.
# The bursts are the dense vertical lines of dots.
burst_timestamps = [2, 11, 17, 26, 36, 47]
print(f"Identified burst timestamps (s): {burst_timestamps}")

# Step 2: Calculate the inter-burst intervals (IBIs).
# IBI is the time difference between consecutive bursts.
inter_burst_intervals = []
for i in range(1, len(burst_timestamps)):
    interval = burst_timestamps[i] - burst_timestamps[i-1]
    inter_burst_intervals.append(interval)

print(f"Calculated inter-burst intervals (s): {inter_burst_intervals}")

# Step 3: Calculate the average inter-burst interval.
average_ibi = np.mean(inter_burst_intervals)

# Display the calculation process
intervals_sum_str = " + ".join(map(str, inter_burst_intervals))
print(f"\nCalculation of the average inter-burst interval:")
print(f"Average = ({intervals_sum_str}) / {len(inter_burst_intervals)}")
print(f"Average = {sum(inter_burst_intervals)} / {len(inter_burst_intervals)}")

# Step 4: Present the final answer to two significant figures.
# Using a format string to ensure two significant figures.
# For a number like 9.0, '{:.2g}'.format(9.0) -> '9', which is not what we want.
# We want to show the precision. '{:.1f}' will show one decimal place.
final_answer = f"{average_ibi:.1f}"
print(f"The average inter-burst interval is {final_answer} s.")
