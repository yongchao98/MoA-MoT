import math

# Step 1: Calculate the total information Elsa can transmit.
video_duration_min = 5
seconds_per_min = 60
blink_rate_per_sec = 1  # 1 blink per second means she can encode 1 bit of info per second

total_seconds = video_duration_min * seconds_per_min
# Each second, she can either blink (1) or not blink (0), which is 1 bit of information.
total_information_bits = total_seconds * blink_rate_per_sec

print(f"The video is {video_duration_min} minutes long, which is {total_seconds} seconds.")
print(f"With a maximum blink rate of {blink_rate_per_sec} per second, Elsa can transmit a total of {total_information_bits} bits.\n")

# Step 2: Calculate the bits required to specify one bomb location.
map_width = 2000
map_height = 1000
total_locations = map_width * map_height

# We need to find 'b' such that 2^b >= total_locations.
# This is calculated as the ceiling of log base 2 of the total locations.
bits_per_location = math.ceil(math.log2(total_locations))

print(f"The map resolution is {map_width}x{map_height}, so there are {total_locations:,} possible locations.")
print(f"To uniquely specify one location, Elsa needs {bits_per_location} bits of information (since 2^{bits_per_location-1} is too small and 2^{bits_per_location} is large enough).\n")

# Step 3: Calculate the maximum number of bomb locations she can express.
# We can only send a whole number of locations, so we use integer division (floor).
max_locations = total_information_bits // bits_per_location

print("To find the maximum number of locations, we divide the total information by the bits required per location.")
print(f"The final calculation is: {total_information_bits} bits / {bits_per_location} bits per location.")
print(f"Result: {max_locations} locations.")
print("\nSo, the highest number of bomb locations Elsa can express is 14.")

<<<14>>>