# 1. Calculate the distance from the repeating pattern over the first 6 minutes.
# There are 4 blocks of 1.5 minutes in the first 6 minutes.
num_blocks = 4
# Each 1.5-minute block yields 2 meters of travel.
dist_per_block = 2
dist_first_6_min = num_blocks * dist_per_block

# 2. Calculate the distance from the final minute.
# The last minute [6, 7] is covered by a single observer, adding 1 meter.
dist_last_1_min = 1

# 3. Calculate the total maximal distance.
total_dist = dist_first_6_min + dist_last_1_min

# 4. Print the final equation.
print(f"The calculation for the maximal distance is:")
print(f"({num_blocks} blocks * {dist_per_block} meters/block) + {dist_last_1_min} meter (for the last minute) = {total_dist} meters")
print(f"{num_blocks} * {dist_per_block} + {dist_last_1_min} = {total_dist}")
