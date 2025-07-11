import math

# --- Parameters ---
video_duration_min = 5
map_width = 2000
map_height = 1000
blink_rate_per_sec = 1

# --- Calculations ---

# Step 1: Calculate the total time in seconds
total_seconds = video_duration_min * 60

# Step 2: Calculate the total bits Elsa can transmit.
# Each second can be a blink (1) or no blink (0), which is 1 bit of information.
total_bits_available = total_seconds * blink_rate_per_sec

# Step 3: Calculate the number of bits needed to specify one bomb location.
# To represent N unique values, we need ceil(log2(N)) bits.
# Bits for x-coordinate (values 0 to 1999, so 2000 possibilities)
bits_for_x = math.ceil(math.log2(map_width))
# Bits for y-coordinate (values 0 to 999, so 1000 possibilities)
bits_for_y = math.ceil(math.log2(map_height))
# Total bits per location is the sum of bits for x and y.
bits_per_location = bits_for_x + bits_for_y

# Step 4: Calculate the maximum number of locations.
# We can only communicate a whole number of locations, so we use floor division.
max_locations = total_bits_available // bits_per_location

# --- Output ---

print("To find the highest number of bomb locations, we follow these steps:\n")

print(f"1. Total information Elsa can transmit:")
print(f"   - Video duration: {video_duration_min} minutes = {video_duration_min} * 60 = {total_seconds} seconds.")
print(f"   - Information rate: {blink_rate_per_sec} blink/second = {blink_rate_per_sec} bit/second.")
print(f"   - Total available bits: {total_seconds} * {blink_rate_per_sec} = {total_bits_available} bits.\n")

print(f"2. Information needed for one location (x, y):")
print(f"   - Bits for x-coordinate (0-{map_width-1}): ceil(log2({map_width})) = {bits_for_x} bits.")
print(f"   - Bits for y-coordinate (0-{map_height-1}): ceil(log2({map_height})) = {bits_for_y} bits.")
print(f"   - Total bits per location: {bits_for_x} + {bits_for_y} = {bits_per_location} bits.\n")

print(f"3. Highest number of locations:")
print(f"   - Calculation: floor(Total bits / Bits per location)")
print(f"   - Equation: floor({total_bits_available} / {bits_per_location}) = {max_locations}\n")

print(f"The highest number of bomb locations Elsa can express is {max_locations}.")
<<<14>>>