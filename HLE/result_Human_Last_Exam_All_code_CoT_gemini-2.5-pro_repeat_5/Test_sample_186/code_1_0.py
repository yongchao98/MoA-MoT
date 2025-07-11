import math

# Define the initial parameters
video_duration_min = 5
fps = 30
blink_rate_per_sec = 1
map_width = 2000
map_height = 1000

# Step 1: Calculate the total communication capacity in bits.
seconds_per_minute = 60
total_seconds = video_duration_min * seconds_per_minute
# Since Elsa can blink or not blink each second, this represents a binary choice.
# Therefore, the total information capacity is equal to the total number of seconds.
total_bits_available = total_seconds

print("Step 1: Calculate the total communication capacity in bits.")
print(f"The video is {video_duration_min} minutes long, which is {video_duration_min} * {seconds_per_minute} = {total_seconds} seconds.")
print("Since Elsa can blink or not blink once per second, she can transmit 1 bit of information per second.")
print(f"Total information capacity: {total_bits_available} bits.\n")

# Step 2: Calculate the bits required to encode one bomb location.
# To represent a number 'n', we need ceil(log2(n)) bits.
bits_for_x = math.ceil(math.log2(map_width))
bits_for_y = math.ceil(math.log2(map_height))
bits_per_location = bits_for_x + bits_for_y

print("Step 2: Calculate the bits required to encode one bomb location.")
print(f"The map resolution is {map_width}x{map_height}.")
print(f"To encode the x-coordinate (from 0 to {map_width - 1}), we need ceil(log2({map_width})) = {bits_for_x} bits.")
print(f"To encode the y-coordinate (from 0 to {map_height - 1}), we need ceil(log2({map_height})) = {bits_for_y} bits.")
print(f"Total bits per location = {bits_for_x} + {bits_for_y} = {bits_per_location} bits.\n")

# Step 3: Calculate the maximum number of bomb locations.
# We use integer division since we can't communicate a fraction of a location.
max_locations = total_bits_available // bits_per_location

print("Step 3: Calculate the maximum number of bomb locations.")
print("Divide the total information capacity by the bits needed per location.")
# Final equation
print(f"Maximum number of locations = floor({total_bits_available} / {bits_per_location}) = {max_locations}.")
print(f"\nTherefore, the highest number of bomb locations Elsa can express is {max_locations}.")
