import math

# Step 1: Calculate total video duration in seconds
video_duration_minutes = 5
seconds_per_minute = 60
total_seconds = video_duration_minutes * seconds_per_minute

# Step 2: Determine total information capacity in bits
# 1 blink per second can be modeled as 1 bit per second
total_bits_capacity = total_seconds

# Step 3: Calculate the information needed for one location
map_resolution_x = 2000
map_resolution_y = 1000

# Bits needed for x-coordinate (0 to 1999)
bits_for_x = math.ceil(math.log2(map_resolution_x))

# Bits needed for y-coordinate (0 to 999)
bits_for_y = math.ceil(math.log2(map_resolution_y))

# Total bits per bomb location
bits_per_location = bits_for_x + bits_for_y

# Step 4: Calculate the maximum number of locations
# We can only communicate a whole number of locations, so we take the floor.
max_locations = math.floor(total_bits_capacity / bits_per_location)

# Print the breakdown of the calculation
print(f"Total video time: {video_duration_minutes} minutes = {total_seconds} seconds")
print(f"Total information Elsa can transmit: {total_bits_capacity} bits")
print("-" * 30)
print(f"Bits to encode x-coordinate (0-{map_resolution_x - 1}): {bits_for_x}")
print(f"Bits to encode y-coordinate (0-{map_resolution_y - 1}): {bits_for_y}")
print(f"Total bits per bomb location: {bits_for_x} + {bits_for_y} = {bits_per_location} bits")
print("-" * 30)
print(f"Maximum number of locations = floor({total_bits_capacity} / {bits_per_location})")
print(f"Final Answer: The highest number of bomb locations Elsa can express is {max_locations}.")
