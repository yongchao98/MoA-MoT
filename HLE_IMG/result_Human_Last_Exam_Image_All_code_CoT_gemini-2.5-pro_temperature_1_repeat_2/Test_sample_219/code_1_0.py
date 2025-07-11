import numpy as np

# Step 1: Identify the timestamps of synchronized bursts from the graph.
# By visual inspection, the dense vertical bands (bursts) are centered at approximately:
burst_timestamps = [2.5, 16.5, 35.5, 47.5]
print(f"Identified burst timestamps (s): {burst_timestamps}")
print("-" * 30)

# Step 2: Calculate the inter-burst intervals (IBIs).
# The IBI is the difference between consecutive burst times.
intervals = np.diff(burst_timestamps)
print("Calculating inter-burst intervals:")
for i in range(len(intervals)):
    print(f"Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {intervals[i]:.1f} s")
print("-" * 30)

# Step 3: Calculate the average inter-burst interval.
average_interval = np.mean(intervals)
print("Calculating the average inter-burst interval:")
interval_sum_str = " + ".join([f"{val:.1f}" for val in intervals])
print(f"Average = ({interval_sum_str}) / {len(intervals)}")
print(f"Average = {np.sum(intervals):.1f} / {len(intervals)} = {average_interval:.1f} s")
print("-" * 30)

# Step 4: Round the result to two significant figures.
# The format specifier 'g' is used for general format, which works well for significant figures.
# In this case, 15.0 becomes 15 (two significant figures).
final_answer = float(f"{average_interval:.2g}")
print(f"The average inter-burst interval rounded to two significant figures is: {final_answer} s")
