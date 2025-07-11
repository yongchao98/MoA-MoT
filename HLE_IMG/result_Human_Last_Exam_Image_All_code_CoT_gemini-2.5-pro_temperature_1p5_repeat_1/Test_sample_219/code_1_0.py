import numpy as np

# Step 1: Identify the timestamps of the synchronized bursts from the graph.
# The bursts are dense vertical lines of dots, indicating many electrodes firing at once.
burst_timestamps = [3, 17, 36, 49]
print(f"Identified burst timestamps (s): {burst_timestamps}")
print("-" * 30)

# Step 2: Calculate the inter-burst intervals (time between consecutive bursts).
intervals = []
print("Calculating inter-burst intervals:")
for i in range(len(burst_timestamps) - 1):
    interval = burst_timestamps[i+1] - burst_timestamps[i]
    intervals.append(interval)
    print(f"Interval {i+1}: {burst_timestamps[i+1]}s - {burst_timestamps[i]}s = {interval}s")
print("-" * 30)

# Step 3: Calculate the average inter-burst interval.
average_interval = np.mean(intervals)
print("Calculating the average inter-burst interval:")
# Building the equation string for clear output
equation_str = f"({'+'.join(map(str, intervals))}) / {len(intervals)}"
print(f"Average = {equation_str} = {average_interval:.4f}s")
print("-" * 30)

# Step 4: Round the average to two significant figures.
# We can use format specifiers for this. '.2g' formats to two significant figures.
avg_rounded = float(f"{average_interval:.2g}")

print(f"The average inter-burst interval rounded to two significant figures is: {avg_rounded}s")