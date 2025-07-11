# The total duration of the snail's journey in minutes.
total_duration = 7

# We use a pattern that gains 2 meters every 1.5 minutes.
# This pattern can be repeated as long as a full 1.5-minute block fits.
block_duration = 1.5
distance_per_block = 2

# The last minute of the journey, [6, 7], must be covered by an observer watching
# that specific interval. This means the snail can only travel 1 meter during it.
final_minute_duration = 1
distance_in_final_minute = 1

# Calculate the duration covered by the repeating blocks.
# This is the total duration minus the final minute.
duration_for_blocks = total_duration - final_minute_duration

# Calculate how many full 1.5-minute blocks can fit into this duration.
num_blocks = int(duration_for_blocks / block_duration)

# Calculate the total distance covered in these blocks.
distance_from_blocks = num_blocks * distance_per_block

# The maximal distance is the sum of the distance from the blocks and
# the distance from the final minute.
max_distance = distance_from_blocks + distance_in_final_minute

# Print the calculation step-by-step
print(f"To find the maximal distance, we can break down the {total_duration}-minute journey:")
print(f"1. A repeating pattern is used for the first {int(duration_for_blocks)} minutes.")
print(f"   - The optimal pattern covers {distance_per_block} meters in {block_duration} minutes.")
print(f"   - Number of blocks that fit: {int(duration_for_blocks)} min / {block_duration} min/block = {num_blocks} blocks.")
print(f"   - Distance from these blocks: {num_blocks} blocks * {distance_per_block} m/block = {distance_from_blocks} meters.")
print(f"\n2. For the final minute ([{int(duration_for_blocks)}, {total_duration}]), the snail's movement is restricted.")
print(f"   - An observer on this interval means the snail advances exactly {distance_in_final_minute} meter.")
print(f"\n3. The total maximal distance is the sum of these parts.")
print(f"   - Final Equation: {distance_from_blocks} meters (from blocks) + {distance_in_final_minute} meter (from last minute) = {max_distance} meters.")
print(f"\nTherefore, the maximal distance the snail could have traveled is {max_distance} meters.")
