import math

# Define the constants based on the problem statement
video_duration_min = 5
map_width = 2000
map_height = 1000
blink_rate_per_sec = 1

# 1. Calculate the total number of bits Elsa can transmit.
# The total time in seconds is the video duration in minutes * 60.
# The data rate is 1 bit/sec (1 blink/sec), so total bits = total seconds.
total_bits = video_duration_min * 60

# 2. Calculate the number of bits required to specify a single bomb location.
# To represent N distinct values, we need ceil(log2(N)) bits.
bits_for_x = math.ceil(math.log2(map_width))
bits_for_y = math.ceil(math.log2(map_height))
bits_per_location = bits_for_x + bits_for_y

# 3. Calculate the maximum number of locations.
# Use integer division because we can only communicate a whole number of locations.
max_locations = total_bits // bits_per_location

# 4. Print the breakdown and the final equation with all numbers.
print(f"Elsa can transmit {total_bits} total bits in {video_duration_min} minutes.")
print(f"Each location requires {bits_for_x} bits for x and {bits_for_y} bits for y, for a total of {bits_per_location} bits.")
print("\nFinal Equation:")
print(f"{total_bits} / ({bits_for_x} + {bits_for_y}) = {max_locations}")