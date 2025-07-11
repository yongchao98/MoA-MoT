import numpy as np

# Step 1: Identify the timestamps of the synchronized bursts from the graph.
# These are estimated from the center of the dense vertical bands of activity.
burst_timestamps = [3.5, 16, 27, 37, 49]
print(f"Identified burst timestamps (s): {burst_timestamps}")

# Step 2: Calculate the inter-burst intervals (IBIs).
# The IBI is the difference in time between consecutive bursts.
inter_burst_intervals = np.diff(burst_timestamps)
print("\nCalculating inter-burst intervals:")
for i in range(len(inter_burst_intervals)):
    print(f"Interval {i+1}: {burst_timestamps[i+1]} - {burst_timestamps[i]} = {inter_burst_intervals[i]} s")

# Step 3: Calculate the average inter-burst interval.
average_ibi = np.mean(inter_burst_intervals)

# Build the equation string for the final output
intervals_sum_str = " + ".join(map(str, inter_burst_intervals))
print("\nCalculating the average inter-burst interval:")
print(f"Average = ({intervals_sum_str}) / {len(inter_burst_intervals)}")
print(f"Average = {sum(inter_burst_intervals)} / {len(inter_burst_intervals)} = {average_ibi:.4f} s")

# Step 4: Round the average to two significant figures.
# We can use a format specifier for this.
rounded_average_ibi = float(f"{average_ibi:.2g}")
print(f"\nRounding the average to two significant figures: {rounded_average_ibi} s")

# The final answer in the required format
# <<<answer content>>>