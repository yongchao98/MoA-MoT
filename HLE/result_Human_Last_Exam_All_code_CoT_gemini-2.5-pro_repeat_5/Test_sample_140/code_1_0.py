import math

def calculate_blocks_in_cuboid(cuboid_dims, block_dims):
    """
    Calculates the maximum number of blocks that can fit in a cuboid.
    It checks all 3 unique orientations for a block like (l, w, w).
    """
    # Sort dimensions to handle permutations easily
    c_l, c_w, c_h = sorted(cuboid_dims, reverse=True)
    b_l, b_w, b_h = sorted(block_dims, reverse=True)

    max_blocks = 0
    
    # Orientation 1: Block (l, w, h) aligns with Cuboid (L, W, H)
    if c_l >= b_l and c_w >= b_w and c_h >= b_h:
        count1 = (c_l // b_l) * (c_w // b_w) * (c_h // b_h)
        max_blocks = max(max_blocks, count1)

    # Orientation 2: Block (l, w, h) aligns with Cuboid (W, L, H) etc.
    # Check all 6 permutations of cuboid dimensions against the block dimensions.
    # Since block is 9x1x1, many permutations are redundant, but this is a general approach.
    import itertools
    for p in set(itertools.permutations(cuboid_dims)):
        c_p_l, c_p_w, c_p_h = p
        if c_p_l >= b_l and c_p_w >= b_w and c_p_h >= b_h:
            count = (c_p_l // b_l) * (c_p_w // b_w) * (c_p_h // b_h)
            max_blocks = max(max_blocks, count)
            
    return max_blocks

# Define dimensions
CUBE_EDGE = 15
BLOCK_DIMS = [9, 1, 1]

# Partition the cube edge
part1 = BLOCK_DIMS[0]
part2 = CUBE_EDGE - part1

# Define the sub-cuboids based on the (9, 6) partition
cuboid_999_dims = (part1, part1, part1)  # 9x9x9
cuboid_996_dims = (part1, part1, part2)  # 9x9x6
cuboid_966_dims = (part1, part2, part2)  # 9x6x6
cuboid_666_dims = (part2, part2, part2)  # 6x6x6

# Calculate how many blocks fit in each type of sub-cuboid
blocks_in_999 = calculate_blocks_in_cuboid(cuboid_999_dims, BLOCK_DIMS)
blocks_in_996 = calculate_blocks_in_cuboid(cuboid_996_dims, BLOCK_DIMS)
blocks_in_966 = calculate_blocks_in_cuboid(cuboid_966_dims, BLOCK_DIMS)
blocks_in_666 = calculate_blocks_in_cuboid(cuboid_666_dims, BLOCK_DIMS)

# There is 1 cuboid of 9x9x9, 3 of 9x9x6, 3 of 9x6x6, and 1 of 6x6x6
total_from_999 = 1 * blocks_in_999
total_from_996 = 3 * blocks_in_996
total_from_966 = 3 * blocks_in_966
total_from_666 = 1 * blocks_in_666

# Calculate the final total
total_blocks = total_from_999 + total_from_996 + total_from_966 + total_from_666

# Print the step-by-step calculation
print(f"The 15x15x15 cube is partitioned into sub-cuboids.")
print(f"Number of blocks in the single 9x9x9 cuboid: {blocks_in_999}")
print(f"Number of blocks in each of the three 9x9x6 cuboids: {blocks_in_996}")
print(f"Number of blocks in each of the three 9x6x6 cuboids: {blocks_in_966}")
print(f"Number of blocks in the single 6x6x6 cuboid: {blocks_in_666}")
print("\nThe final calculation is:")
print(f"(1 * {blocks_in_999}) + (3 * {blocks_in_996}) + (3 * {blocks_in_966}) + (1 * {blocks_in_666})")
print(f"= {total_from_999} + {total_from_996} + {total_from_966} + {total_from_666}")
print(f"= {total_blocks}")
print(f"\nThe largest number of blocks that will fit is {total_blocks}.")
