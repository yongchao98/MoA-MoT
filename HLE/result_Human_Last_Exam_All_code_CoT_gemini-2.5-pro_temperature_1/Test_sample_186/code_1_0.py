import math

# Step 1: Calculate total communication time in seconds
video_duration_min = 5
seconds_per_minute = 60
total_seconds = video_duration_min * seconds_per_minute

print(f"Video duration is {video_duration_min} minutes, which is {total_seconds} seconds.")

# Step 2: Determine the total information bandwidth in bits
# Elsa can blink once per second, which can be used to encode 1 bit of data per second.
# A blink = '1', no blink = '0'.
bits_per_second = 1
total_bits_available = total_seconds * bits_per_second

print(f"With a rate of {bits_per_second} blink/second, Elsa can transmit {total_bits_available} bits in total.")
print("-" * 20)

# Step 3: Calculate the number of bits needed to encode one bomb location
map_width = 2000  # values from 0 to 1999
map_height = 1000 # values from 0 to 999

# We need to find the number of bits to represent all possible values for each coordinate.
# The formula is ceil(log2(N)), where N is the number of values.
bits_for_x = math.ceil(math.log2(map_width))
bits_for_y = math.ceil(math.log2(map_height))
bits_per_location = bits_for_x + bits_for_y

print(f"To represent a width of {map_width}, we need {bits_for_x} bits.")
print(f"To represent a height of {map_height}, we need {bits_for_y} bits.")
print(f"Therefore, to express one bomb location ({bits_for_x} + {bits_for_y}), it takes {bits_per_location} bits (or seconds).")
print("-" * 20)

# Step 4: Calculate the highest number of bomb locations
# We can't send a fraction of a location, so we use integer division.
max_locations = total_bits_available // bits_per_location

print("Final Calculation:")
print(f"The total number of locations is the total available bits divided by the bits per location.")
print(f"{total_bits_available} (total bits) // {bits_per_location} (bits per location) = {max_locations}")
print("-" * 20)
print(f"The highest number of bomb locations Elsa can express is {max_locations}.")
