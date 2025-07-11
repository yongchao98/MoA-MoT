import numpy as np

# 1. Identify the timestamps of the synchronized bursts from the graph.
# These are the approximate time points (in seconds) where dense vertical firing patterns occur.
burst_timestamps = [3, 15, 24, 36, 49]
print(f"Identified burst timestamps (s): {burst_timestamps}\n")

# 2. Calculate the inter-burst intervals (IBIs).
# The interval is the time difference between consecutive bursts.
inter_burst_intervals = np.diff(burst_timestamps)
print(f"Calculated inter-burst intervals (s): {list(inter_burst_intervals)}\n")

# 3. Calculate the average inter-burst interval.
average_interval = np.mean(inter_burst_intervals)
sum_of_intervals = np.sum(inter_burst_intervals)
num_of_intervals = len(inter_burst_intervals)

# 4. Print the calculation steps and the final result.
# Create a string representation of the sum for the equation.
interval_sum_str = " + ".join(map(str, inter_burst_intervals))

print("Calculation of the average inter-burst interval:")
# Print the full equation for calculating the average
print(f"({interval_sum_str}) / {num_of_intervals} = {sum_of_intervals} / {num_of_intervals} = {average_interval} s\n")

# Round the average to two significant figures.
# 11.5 rounds to 12 when rounding to two significant figures.
rounded_average = float(f"{average_interval:.2g}")
print(f"The average inter-burst interval rounded to two significant figures is: {rounded_average} s")