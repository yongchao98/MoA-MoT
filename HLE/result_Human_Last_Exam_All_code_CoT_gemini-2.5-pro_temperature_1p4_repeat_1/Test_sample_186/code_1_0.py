import math

# Define the given parameters
video_duration_min = 5
map_width = 2000
map_height = 1000

# 1. Calculate the total time in seconds, which corresponds to the total bits of information.
total_seconds = video_duration_min * 60
total_bits = total_seconds  # 1 blink per second = 1 bit per second

# 2. Calculate the number of bits needed to encode a single location (x, y).
# Bits for x-coordinate (values 0 to 1999)
bits_for_x = math.ceil(math.log2(map_width))
# Bits for y-coordinate (values 0 to 999)
bits_for_y = math.ceil(math.log2(map_height))
# Total bits for one location
bits_per_location = bits_for_x + bits_for_y

# 3. Calculate the maximum number of complete locations that can be communicated.
max_locations = math.floor(total_bits / bits_per_location)

# --- Output the result with the final equation ---
print("This problem can be solved by calculating the total information Elsa can transmit and dividing it by the information needed for a single bomb location.")
print("\nFirst, we calculate the total bits Elsa can transmit:")
print(f"Total bits = {video_duration_min} min * 60 sec/min * 1 bit/sec = {total_bits} bits")

print("\nNext, we calculate the bits needed for one location's coordinates:")
print(f"Bits for x (0-{map_width-1}): ceil(log2({map_width})) = {bits_for_x} bits")
print(f"Bits for y (0-{map_height-1}): ceil(log2({map_height})) = {bits_for_y} bits")
print(f"Total bits per location = {bits_for_x} + {bits_for_y} = {bits_per_location} bits")

print("\nFinally, we can find the highest number of locations:")
print("Highest number of locations = floor(Total bits / Bits per location)")
print(f"Result = floor({total_bits} / {bits_per_location}) = {max_locations}")
<<<14>>>