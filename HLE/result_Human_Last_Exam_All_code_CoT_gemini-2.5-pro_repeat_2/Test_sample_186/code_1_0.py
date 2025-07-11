import math

# Step 1: Calculate Total Information Capacity
video_duration_minutes = 5
seconds_per_minute = 60
blink_rate_per_second = 1 # This is our bit rate

# Total time in seconds
total_seconds = video_duration_minutes * seconds_per_minute

# Total bits Elsa can transmit (1 bit per second)
total_bits_capacity = total_seconds * blink_rate_per_second

# Step 2: Calculate Information per Location
map_width = 2000
map_height = 1000

# Total number of unique locations on the map
total_possible_locations = map_width * map_height

# Bits needed to specify one unique location
# This is calculated as ceil(log2(N)) where N is the number of locations
bits_per_location = math.ceil(math.log2(total_possible_locations))

# Step 3: Calculate Maximum Number of Locations
# We use integer division because Elsa cannot send information for a partial location
max_locations = total_bits_capacity // bits_per_location

# Output the explanation and the final calculation
print(f"Total communication time: {video_duration_minutes} min * {seconds_per_minute} s/min = {total_seconds} seconds.")
print(f"Information can be sent at {blink_rate_per_second} bit/second, so total capacity is {total_bits_capacity} bits.")
print(f"Total possible bomb locations on a {map_width}x{map_height} map: {total_possible_locations}.")
print(f"Bits needed to encode one location: ceil(log2({total_possible_locations})) = {bits_per_location} bits.")
print("\nTo find the maximum number of locations, we divide the total capacity by the bits per location:")
print(f"Final Equation: {total_bits_capacity} // {bits_per_location} = {max_locations}")
print(f"\nTherefore, the highest number of bomb locations Elsa can express is {max_locations}.")
