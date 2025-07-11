import math

# Step 1: Define the problem's parameters
video_duration_minutes = 5
blinks_per_second = 1
map_width = 2000
map_height = 1000

# Step 2: Calculate the total number of seconds in the video.
# This represents the total number of bits Elsa can transmit.
total_seconds = video_duration_minutes * 60
total_bits_available = total_seconds * blinks_per_second

print(f"Elsa can transmit for a total of {total_seconds} seconds.")
print(f"This allows her to send {total_bits_available} bits of information.")
print("-" * 20)

# Step 3: Calculate the total number of unique locations on the map.
total_possible_locations = map_width * map_height

print(f"The map has {total_possible_locations:,} unique locations.")

# Step 4: Calculate the number of bits required to specify a single location.
# We need to find the smallest integer 'b' such that 2^b >= total_possible_locations.
bits_per_location = math.ceil(math.log2(total_possible_locations))

print(f"Each location requires {bits_per_location} bits to be uniquely identified.")
print("-" * 20)

# Step 5: Calculate the highest number of bomb locations Elsa can express.
# This is the total bits available divided by the bits needed per location.
# We use integer division because she can only communicate a whole number of locations.
max_locations = total_bits_available // bits_per_location

print("To find the maximum number of locations, we perform the following calculation:")
print(f"Maximum Locations = Total Bits Available / Bits Per Location")
print(f"Maximum Locations = {total_bits_available} / {bits_per_location}")
print(f"Result: {max_locations}")

# The final answer in the required format
# <<<14>>>