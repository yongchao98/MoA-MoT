# Plan:
# The maximum rate of travel the snail can achieve is 2 meters per 1.5 minutes.
# We apply this rate for as many full 1.5-minute segments as possible in the 7-minute total time.
# Any remaining time is covered by a standard 1-meter-per-minute rate.

total_time = 7  # minutes
segment_time = 1.5  # minutes per optimal segment
segment_distance = 2  # meters per optimal segment

# Calculate how many full optimal segments fit into the total time.
num_segments = int(total_time / segment_time)

# Calculate the time covered by these optimal segments.
time_in_segments = num_segments * segment_time

# Calculate the distance covered in these segments.
distance_from_segments = num_segments * segment_distance

# Calculate the remaining time.
remaining_time = total_time - time_in_segments

# The distance covered in the remaining time is 1 meter per minute.
# Since each observer watches for exactly 1 minute and sees 1 meter traveled,
# the distance covered in the remaining_time duration is equal to the duration itself.
distance_from_remainder = remaining_time

# The total maximum distance is the sum of the distances from all parts.
total_distance = distance_from_segments + distance_from_remainder

print("The calculation for the maximal distance is as follows:")
print(f"Number of {segment_time}-minute segments: {num_segments}")
print(f"Distance from these segments: {num_segments} * {segment_distance} = {distance_from_segments} meters")
print(f"Remaining time: {total_time} - ({num_segments} * {segment_time}) = {remaining_time} minute(s)")
print(f"Distance from remaining time: {distance_from_remainder} meter(s)")
print(f"Final equation: {distance_from_segments} + {distance_from_remainder} = {total_distance}")
print(f"\nThe maximal distance the snail could have advanced is {total_distance} meters.")
