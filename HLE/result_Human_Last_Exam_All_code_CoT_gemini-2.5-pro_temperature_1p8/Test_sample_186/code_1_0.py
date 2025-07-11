import math

# --- Problem Parameters ---
video_duration_minutes = 5
map_width = 2000
map_height = 1000
blinks_per_second = 1 # This implies 1 bit of information can be sent per second

# --- Step 1: Calculate total communication time in seconds ---
seconds_per_minute = 60
total_seconds = video_duration_minutes * seconds_per_minute

# --- Step 2: Determine total information capacity in bits ---
# A blink or no-blink in a 1-second interval can represent 1 bit of data.
total_bits_available = total_seconds * blinks_per_second

# --- Step 3: Calculate the information required per location ---
# First, find the total number of possible locations
total_locations = map_width * map_height
# Then, find the bits needed to represent one location.
# This is ceil(log2(total_locations)) because we need to cover all possibilities
# and we can't use a fraction of a bit.
bits_per_location = math.ceil(math.log2(total_locations))

# --- Step 4: Calculate the maximum number of locations ---
# This is the total information we can send divided by the information per message.
# We use integer division because she cannot communicate a fraction of a location.
max_bomb_locations = total_bits_available // bits_per_location

# --- Final Output ---
print("This script calculates the maximum number of bomb locations Elsa can communicate.")
print("---")
print(f"1. Total communication time: {video_duration_minutes} min * {seconds_per_minute} s/min = {total_seconds} seconds.")
print(f"2. Total information capacity: {total_seconds} seconds * {blinks_per_second} bit/sec = {total_bits_available} bits.")
print(f"3. Information needed per location: ceil(log2({map_width} * {map_height})) = {bits_per_location} bits.")
print("---")
print("Final Calculation:")
print(f"To find the maximum number of locations, we perform the calculation: {total_bits_available} (total bits) / {bits_per_location} (bits per location)")
print("The final equation for the maximum number of whole locations is:")
print(f"{max_bomb_locations} = {total_bits_available} // {bits_per_location}")
print(f"\nThe highest number of bomb locations can Elsa expresses is: {max_bomb_locations}")