import numpy as np

# 1. Identify the timestamps of the synchronized bursts from the graph.
burst_timestamps = [3, 17, 24, 36, 49]
print(f"Identified burst timestamps (s): {burst_timestamps}")

# 2. Calculate the inter-burst intervals (IBIs).
inter_burst_intervals = []
for i in range(len(burst_timestamps) - 1):
  interval = burst_timestamps[i+1] - burst_timestamps[i]
  inter_burst_intervals.append(interval)

print(f"Calculated inter-burst intervals (s): {inter_burst_intervals}")

# 3. Calculate the average inter-burst interval.
average_ibi = np.mean(inter_burst_intervals)

# 4. Print the final calculation and result.
# Create the sum part of the equation string
sum_str = " + ".join(map(str, inter_burst_intervals))
# Get the number of intervals
num_intervals = len(inter_burst_intervals)

print("\nCalculating the average inter-burst interval:")
print(f"Average = ({sum_str}) / {num_intervals}")
print(f"Average = {sum(inter_burst_intervals)} / {num_intervals}")
print(f"Average = {average_ibi:.3f} s")

# 5. Round the result to two significant figures.
# A simple way to round to 2 significant figures for a number >= 10 is to round to the nearest integer.
# For a number < 10, it would be rounding to one decimal place.
# Since 11.5 >= 10, we round to zero decimal places.
final_answer = round(average_ibi)
print(f"\nThe average inter-burst interval rounded to two significant figures is {final_answer} s.")
