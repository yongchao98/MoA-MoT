import math

# Step 1: Define the constants from the problem description.
video_duration_minutes = 5
seconds_per_minute = 60
blinks_per_second = 1
map_width = 2000
map_height = 1000

# Step 2: Calculate the total information capacity in bits.
# At each second, Elsa can either blink (1) or not blink (0), which represents one bit of information.
total_seconds = video_duration_minutes * seconds_per_minute
total_bits_available = total_seconds * blinks_per_second

print("### Step 1: Calculate Total Information Capacity ###")
print(f"Total time available: {video_duration_minutes} minutes * {seconds_per_minute} seconds/minute = {total_seconds} seconds")
print(f"Total information capacity: {total_seconds} seconds * {blinks_per_second} bit/second = {total_bits_available} bits")
print("-" * 30)

# Step 3: Calculate the number of bits required to encode one bomb location (x, y).
# We use ceiling of log base 2 to find the bits needed for each coordinate.
bits_for_x = math.ceil(math.log2(map_width))
bits_for_y = math.ceil(math.log2(map_height))
bits_per_location = bits_for_x + bits_for_y

print("### Step 2: Calculate Bits Needed Per Location ###")
print(f"Map resolution: {map_width}x{map_height}")
print(f"Bits needed for x-coordinate (up to {map_width}): ceil(log2({map_width})) = {bits_for_x} bits")
print(f"Bits needed for y-coordinate (up to {map_height}): ceil(log2({map_height})) = {bits_for_y} bits")
print(f"Total bits per location: {bits_for_x} + {bits_for_y} = {bits_per_location} bits")
print("-" * 30)

# Step 4: Calculate the maximum number of bomb locations.
# This is the total bits available divided by the bits needed for each location.
# We must take the floor because Elsa can only communicate complete locations.
max_locations = math.floor(total_bits_available / bits_per_location)

print("### Step 3: Calculate Max Number of Bomb Locations ###")
print("This is the total capacity divided by the cost per location.")
print(f"Max locations = floor(Total Bits Available / Bits Per Location)")
print(f"Max locations = floor({total_bits_available} / {bits_per_location}) = {max_locations}")
print("-" * 30)

print(f"\nThe highest number of bomb locations Elsa can express is {max_locations}.")

<<<14>>>