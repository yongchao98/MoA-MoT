import numpy as np

# Step 1: Identify the timestamps of synchronized bursts from the graph.
# Visually inspecting the plot, the dense vertical bands (bursts) are centered at approximately:
burst_timestamps = [2.5, 17.5, 35.5, 48.5]
print(f"Identified burst timestamps (s): {burst_timestamps}")
print("-" * 30)

# Step 2: Calculate the inter-burst intervals (IBIs).
# The IBI is the difference in time between consecutive bursts.
intervals = []
print("Calculating inter-burst intervals:")
for i in range(len(burst_timestamps) - 1):
    interval = burst_timestamps[i+1] - burst_timestamps[i]
    intervals.append(interval)
    print(f"Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {interval:.1f} s")
print("-" * 30)

# Step 3: Calculate the average IBI.
average_ibi = np.mean(intervals)
print("Calculating the average inter-burst interval:")
# Construct the equation string
equation_str = f"({intervals[0]:.1f} + {intervals[1]:.1f} + {intervals[2]:.1f}) / 3 = {average_ibi:.3f} s"
print(equation_str)
print("-" * 30)

# Step 4: Round the result to two significant figures.
# To round to 2 significant figures, we can use a format specifier.
# For a number like 15.333..., two significant figures would be 15.
rounded_average_ibi = float(f"{average_ibi:.2g}")
print(f"The average inter-burst interval rounded to two significant figures is: {rounded_average_ibi} s")
