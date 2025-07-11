import numpy as np

# Step 1 & 2: Identify the timestamps of synchronized bursts from the graph.
# The bursts are dense vertical lines of activity.
burst_timestamps = [3, 16, 26, 36, 48]
print(f"Identified burst timestamps (s): {burst_timestamps}")
print("-" * 30)

# Step 3: Calculate the inter-burst intervals (IBIs).
inter_burst_intervals = []
print("Calculating inter-burst intervals:")
for i in range(1, len(burst_timestamps)):
    interval = burst_timestamps[i] - burst_timestamps[i-1]
    inter_burst_intervals.append(interval)
    print(f"Interval {i}: {burst_timestamps[i]} s - {burst_timestamps[i-1]} s = {interval} s")

print(f"\nCalculated inter-burst intervals (s): {inter_burst_intervals}")
print("-" * 30)

# Step 4: Calculate the average IBI.
average_ibi = np.mean(inter_burst_intervals)
# Build the string for the equation
sum_str = " + ".join(map(str, inter_burst_intervals))

print("Calculating the average inter-burst interval:")
print(f"Average = ({sum_str}) / {len(inter_burst_intervals)}")
print(f"Average = {sum(inter_burst_intervals)} / {len(inter_burst_intervals)} = {average_ibi:.4f} s")
print("-" * 30)

# Step 5: Round the average to two significant figures.
# The format '{:.2g}' handles rounding to the specified number of significant figures.
formatted_average_ibi = "{:.2g}".format(average_ibi)

print(f"The average inter-burst interval is {average_ibi:.4f} s.")
print(f"Rounding to two significant figures, the result is: {formatted_average_ibi} s")

# Final answer in the required format
final_answer = float(formatted_average_ibi)