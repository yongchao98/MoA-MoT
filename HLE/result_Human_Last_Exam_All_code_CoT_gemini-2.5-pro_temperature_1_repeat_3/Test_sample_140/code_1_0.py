def solve_packing_problem():
    """
    Calculates the maximum number of 9x1x1 blocks that fit in a 15x15x15 cube.
    """
    # Define cube and block dimensions
    cube_edge = 15
    block_long = 9
    block_short = 1

    # This problem is a 3D packing problem. We use a greedy partitioning strategy.

    # Step 1: Fill a 9x15x15 section of the 15x15x15 cube.
    # The number of blocks is calculated by dividing each dimension of the section
    # by the corresponding dimension of the block.
    # We align the block's 9-unit length with the 9-unit edge of this section.
    blocks_part1 = (9 // block_long) * (15 // block_short) * (15 // block_short)
    # The remaining space is a (15-9)x15x15 = 6x15x15 cuboid.

    # Step 2: Fill the remaining 6x15x15 space.
    # We must align the block's 9-unit length with a 15-unit edge.
    # We fill a 6x9x15 section of the remaining space.
    blocks_part2 = (6 // block_short) * (9 // block_long) * (15 // block_short)
    # The remaining space is now 6x(15-9)x15 = 6x6x15.

    # Step 3: Fill the remaining 6x6x15 space.
    # We align the block's 9-unit length with the 15-unit edge.
    # We fill a 6x6x9 section.
    blocks_part3 = (6 // block_short) * (6 // block_short) * (9 // block_long)
    # The final remaining space is 6x6x(15-9) = 6x6x6, where no more blocks fit.

    # Calculate the total number of blocks
    total_blocks = blocks_part1 + blocks_part2 + blocks_part3

    # Print the final equation showing how the total is calculated
    print(f"The largest number of blocks is the sum from three partitions.")
    print(f"The calculation is: {blocks_part1} + {blocks_part2} + {blocks_part3} = {total_blocks}")

solve_packing_problem()