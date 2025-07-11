import numpy as np

# Step 1 & 2: Identify the timestamps of the synchronized bursts from the graph.
# The bursts are the dense vertical collections of firing events.
# By inspecting the x-axis (Time (s)), we can estimate their start times.
burst_timestamps = [2, 17, 36, 49]
print(f"Identified burst timestamps (s): {burst_timestamps}")

# Step 3: Calculate the inter-burst intervals (IBIs).
# IBI is the time difference between consecutive bursts.
intervals = np.diff(burst_timestamps)
print(f"Calculated inter-burst intervals (s): {intervals.tolist()}")

# Step 4: Calculate the average IBI.
average_interval = np.mean(intervals)

# Step 5: Display the full calculation and the final answer.
# The format 'g' is used for general number formatting, and '.2' specifies two significant digits.
print(f"The average inter-burst interval is calculated as: ({intervals[0]} + {intervals[1]} + {intervals[2]}) / 3 = {average_interval:.4f} s")
print(f"Rounding to two significant figures, the average inter-burst interval is: {average_interval:.2g} s")