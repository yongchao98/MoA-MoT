import numpy as np

# Step 1: Define the timestamps of the synchronized bursts identified from the graph.
# These are the approximate centers of the dense vertical firing bands.
burst_timestamps = [2.5, 16.5, 34.5, 47.5]
print(f"Identified synchronized burst timestamps: {burst_timestamps} s")

# Step 2: Calculate the inter-burst intervals (IBIs), which are the differences
# between consecutive timestamps.
intervals = np.diff(burst_timestamps)
print("\nCalculating the inter-burst intervals (IBIs):")
for i, interval in enumerate(intervals):
    print(f"Interval {i+1}: {burst_timestamps[i+1]:.1f} s - {burst_timestamps[i]:.1f} s = {interval:.1f} s")

# Step 3: Calculate the average of these intervals.
average_interval = np.mean(intervals)

# Step 4: Display the final calculation and the result rounded to two significant figures.
print("\nCalculating the average inter-burst interval:")
interval_sum_str = " + ".join([f"{i:.1f}" for i in intervals])
print(f"Average IBI = ({interval_sum_str}) s / {len(intervals)}")
print(f"Average IBI = {np.sum(intervals):.1f} s / {len(intervals)}")
print(f"Average IBI = {average_interval:.3f} s")

# The result 15.0 has three significant figures. To two significant figures, it is 15.
# We use the 'g' format specifier to handle this correctly.
final_answer = f"{average_interval:.2g}"
print(f"\nThe average inter-burst interval to two significant figures is: {final_answer} s")
