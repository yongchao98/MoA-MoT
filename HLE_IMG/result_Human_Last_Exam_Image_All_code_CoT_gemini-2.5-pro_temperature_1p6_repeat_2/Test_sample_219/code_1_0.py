import numpy as np

# 1. Identify the timestamps of synchronized bursts from the graph.
# By visual inspection, the center of each major burst occurs at approximately these times in seconds.
burst_timestamps = [4, 15, 25, 36, 47]
print(f"Identified burst timestamps (s): {burst_timestamps}")
print("-" * 30)

# 2. Calculate the inter-burst intervals (IBIs).
# The interval is the time difference between consecutive bursts.
intervals = np.diff(burst_timestamps)
print("Calculating the inter-burst intervals (s):")
for i in range(len(intervals)):
    print(f"Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {intervals[i]} s")
print(f"\nAll intervals: {intervals.tolist()} s")
print("-" * 30)

# 3. Compute the average IBI.
average_ibi = np.mean(intervals)

# 4. Display the final calculation equation.
# We build a string to show how the average was calculated.
numerator_str = " + ".join(map(str, intervals))
equation = f"Average IBI = ({numerator_str}) / {len(intervals)}"
print("Calculation for the average inter-burst interval:")
print(f"{equation} = {average_ibi:.4f} s")
print("-" * 30)

# 5. Round the result to two significant figures.
# For 10.75, the two significant figures are 1 and 0. The next digit is 7,
# so we round up, resulting in 11.
average_ibi_rounded = float(f"{average_ibi:.2g}")
print(f"The final average inter-burst interval, rounded to two significant figures, is {average_ibi_rounded} s.")
