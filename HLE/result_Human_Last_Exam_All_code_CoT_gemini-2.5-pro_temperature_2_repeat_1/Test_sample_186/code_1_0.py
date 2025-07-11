import math

# Step 1: Calculate total video duration in seconds.
video_duration_minutes = 5
seconds_per_minute = 60
total_seconds = video_duration_minutes * seconds_per_minute

# Step 2: Calculate total bits Elsa can transmit.
# A blink or no-blink per second can encode 1 bit of information.
# So, the information rate is 1 bit per second.
total_bits_available = total_seconds * 1

# Step 3: Calculate the number of bits needed per bomb location.
map_width = 2000
map_height = 1000

# Bits for x-coordinate (0 to 1999)
bits_for_x = math.ceil(math.log2(map_width))

# Bits for y-coordinate (0 to 999)
bits_for_y = math.ceil(math.log2(map_height))

# Total bits per location
bits_per_location = bits_for_x + bits_for_y

# Step 4: Calculate the maximum number of locations.
# This is the total available bits divided by the bits needed for one location.
# We take the floor because we can't communicate a fraction of a location.
max_locations = math.floor(total_bits_available / bits_per_location)

print("--- Calculation Steps ---")
print(f"Total video time in seconds: {video_duration_minutes} min * {seconds_per_minute} s/min = {total_seconds} seconds")
print(f"Total information Elsa can transmit (at 1 bit/sec): {total_seconds} bits")
print(f"Bits needed for X-coordinate (up to {map_width}): ceil(log2({map_width})) = {bits_for_x} bits")
print(f"Bits needed for Y-coordinate (up to {map_height}): ceil(log2({map_height})) = {bits_for_y} bits")
print(f"Total bits per bomb location: {bits_for_x} + {bits_for_y} = {bits_per_location} bits")
print("\n--- Final Calculation ---")
print(f"Maximum number of locations = floor(Total bits / Bits per location)")
print(f"Equation: floor({total_bits_available} / {bits_per_location}) = {max_locations}")
print("\nTherefore, the highest number of bomb locations Elsa can express is calculated above.")

print(f"\n<<<{max_locations}>>>")