import numpy as np

# 1. Identify the timestamps of synchronized bursts from the graph.
# By visual inspection, the bursts are centered at these approximate times in seconds.
burst_timestamps = [2.5, 12.0, 18.0, 25.0, 36.0, 47.0]
print(f"Identified burst timestamps (s): {burst_timestamps}")
print("-" * 30)

# 2. Calculate the inter-burst intervals (the time between consecutive bursts).
# This is done by finding the difference between each timestamp and the one before it.
inter_burst_intervals = np.diff(burst_timestamps)
print("Calculated inter-burst intervals (s):")
for i in range(len(inter_burst_intervals)):
    print(f"  Interval {i+1}: {burst_timestamps[i+1]} - {burst_timestamps[i]} = {inter_burst_intervals[i]:.1f} s")
print("-" * 30)

# 3. Calculate the average of these intervals.
average_interval = np.mean(inter_burst_intervals)
sum_of_intervals = np.sum(inter_burst_intervals)
number_of_intervals = len(inter_burst_intervals)

# 4. Print the final calculation and the result rounded to two significant figures.
equation_numerator = ' + '.join([f"{i:.1f}" for i in inter_burst_intervals])
print("Final Calculation:")
print(f"Average Interval = ({equation_numerator}) / {number_of_intervals}")
print(f"Average Interval = {sum_of_intervals:.1f} / {number_of_intervals}")

# The '{:.2g}' format specifier rounds a number to the specified number of significant figures.
print(f"\nThe average inter-burst interval is {average_interval:.2g} s.")
