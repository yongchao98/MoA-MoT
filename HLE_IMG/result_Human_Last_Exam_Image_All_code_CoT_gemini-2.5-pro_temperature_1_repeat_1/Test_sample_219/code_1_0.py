import math

def round_to_sf(num, sf):
    """
    Rounds a number to a specified number of significant figures.
    """
    if num == 0:
        return 0
    # Use 'g' format specifier for significant figures
    return float(f"{num:.{sf}g}")

# 1. By visual inspection of the graph, the synchronized bursts
# (dense vertical bands of activity) are centered at approximately
# 3 s, 17 s, 35 s, and 48 s.
burst_timestamps = [3, 17, 35, 48]
print(f"Identified burst timestamps: {burst_timestamps} s")

# 2. Calculate the time intervals between consecutive bursts.
intervals = []
print("\nCalculating inter-burst intervals (IBIs):")
for i in range(len(burst_timestamps) - 1):
    interval = burst_timestamps[i+1] - burst_timestamps[i]
    intervals.append(interval)
    print(f"Interval {i+1}: {burst_timestamps[i+1]} - {burst_timestamps[i]} = {interval} s")

# 3. Calculate the average of the inter-burst intervals.
sum_of_intervals = sum(intervals)
number_of_intervals = len(intervals)
average_ibi = sum_of_intervals / number_of_intervals

print("\nCalculating the average IBI:")
# Build the equation string with each number
equation_str_numerator = " + ".join(map(str, intervals))
print(f"Average IBI = ({equation_str_numerator}) / {number_of_intervals}")
print(f"Average IBI = {sum_of_intervals} / {number_of_intervals}")
print(f"Average IBI = {average_ibi:.1f} s")

# 4. Round the result to two significant figures.
average_ibi_rounded = round_to_sf(average_ibi, 2)
print(f"\nThe final average inter-burst interval rounded to two significant figures is {average_ibi_rounded} s.")