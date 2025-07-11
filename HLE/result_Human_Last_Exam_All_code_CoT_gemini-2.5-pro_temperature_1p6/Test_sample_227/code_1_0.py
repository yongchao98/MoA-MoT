# Plan:
# 1. Decompose the 7-minute interval into optimized blocks.
#    A 1.5-minute block can yield a maximum distance of 2 meters.
#    A 1-minute block can yield a maximum distance of 1 meter.
# 2. We decompose 7 minutes into four 1.5-minute blocks and one 1-minute block.
#    7 = 1.5 + 1.5 + 1.5 + 1.5 + 1
# 3. Calculate the total distance by summing the contributions from each block.

# Distance from the four 1.5-minute blocks
dist_1_5_min_block = 2
num_1_5_min_blocks = 4
total_dist_from_1_5_min_blocks = dist_1_5_min_block * num_1_5_min_blocks

# Distance from the remaining 1-minute block
dist_1_min_block = 1
num_1_min_blocks = 1
total_dist_from_1_min_blocks = dist_1_min_block * num_1_min_blocks

# The total maximal distance is the sum of the parts.
total_distance = total_dist_from_1_5_min_blocks + total_dist_from_1_min_blocks

# Print the final equation as requested.
distances = [str(dist_1_5_min_block)] * num_1_5_min_blocks + [str(dist_1_min_block)] * num_1_min_blocks
equation = " + ".join(distances)
print(f"The maximal distance is the sum of the distances from each block:")
print(f"{equation} = {total_distance}")