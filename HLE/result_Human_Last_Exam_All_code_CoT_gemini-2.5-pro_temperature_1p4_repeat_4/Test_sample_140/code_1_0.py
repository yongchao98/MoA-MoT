# Plan: Partition the 15x15x15 cube and calculate the number of blocks in each partition.

# Partition 1: 9x9x9 cube
# We align the blocks' 9-unit length with one of the cube's 9-unit dimensions.
# The other dimensions (1x1) fit into the remaining 9x9 area.
blocks_in_9x9x9 = (9 // 1) * (9 // 1) * (9 // 9)
print(f"Blocks in the 9x9x9 region: (9/1) * (9/1) * (9/9) = {blocks_in_9x9x9}")

# Partition 2: Three 9x9x6 slabs
# The 9-unit length of the blocks must align with a 9-unit dimension of the slab.
blocks_in_one_9x9x6 = (9 // 1) * (6 // 1) * (9 // 9)
print(f"Blocks in one 9x9x6 region: (9/1) * (6/1) * (9/9) = {blocks_in_one_9x9x6}")

# Partition 3: Three 9x6x6 rods
# The 9-unit length of the blocks must align with the 9-unit dimension of the rod.
blocks_in_one_9x6x6 = (6 // 1) * (6 // 1) * (9 // 9)
print(f"Blocks in one 9x6x6 region: (6/1) * (6/1) * (9/9) = {blocks_in_one_9x6x6}")

# Partition 4: One 6x6x6 cube
# No dimension is >= 9, so no blocks can fit.
blocks_in_6x6x6 = 0
print(f"Blocks in the 6x6x6 region: {blocks_in_6x6x6}")

# Total calculation
total_blocks = blocks_in_9x9x9 + 3 * blocks_in_one_9x9x6 + 3 * blocks_in_one_9x6x6 + blocks_in_6x6x6

print("\nTotal Calculation:")
print(f"Total = {blocks_in_9x9x9} + 3 * {blocks_in_one_9x9x6} + 3 * {blocks_in_one_9x6x6} + {blocks_in_6x6x6} = {total_blocks}")
