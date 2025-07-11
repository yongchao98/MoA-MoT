import numpy as np

# Step 1: Identify the timestamps of the synchronized bursts from the graph.
# These are the approximate time points (in seconds) at the center of each burst.
burst_timestamps = [3, 16, 26, 36, 49]
print(f"Identified burst timestamps (s): {burst_timestamps}")

# Step 2: Calculate the intervals between consecutive bursts.
inter_burst_intervals = np.diff(burst_timestamps)
print(f"Calculated inter-burst intervals (s): {list(inter_burst_intervals)}")

# Step 3: Calculate the average of these intervals.
average_ibi = np.mean(inter_burst_intervals)

# Step 4: Display the full calculation.
sum_string = " + ".join(map(str, inter_burst_intervals))
print(f"\nCalculating the average inter-burst interval:")
print(f"Average = ({sum_string}) / {len(inter_burst_intervals)}")
print(f"Average = {np.sum(inter_burst_intervals)} / {len(inter_burst_intervals)}")
print(f"Average = {average_ibi:.1f} s")

# Step 5: Round the result to two significant figures.
# 11.5 rounded to two significant figures is 12.
final_answer = 12
print(f"\nRounding to two significant figures, the average inter-burst interval is {final_answer} s.")