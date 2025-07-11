import math

def solve_packing():
    """
    Calculates the maximum number of 9x1x1 blocks that can fit in a 15x15x15 cube.
    """

    # Step 1: Pack a 15x15x9 slab within the 15x15x15 cube.
    # Orient blocks as 1x1x9.
    # Number of blocks = floor(15/1) * floor(15/1) * floor(9/9)
    blocks_1 = math.floor(15/1) * math.floor(15/1) * math.floor(9/9)

    # This leaves a 15x15x6 remaining volume.

    # Step 2: Pack a 15x9x6 slab within the 15x15x6 remaining volume.
    # Orient blocks as 1x9x1.
    # Number of blocks = floor(15/1) * floor(9/9) * floor(6/1)
    blocks_2 = math.floor(15/1) * math.floor(9/9) * math.floor(6/1)

    # This leaves a 15x6x6 remaining volume.

    # Step 3: Pack a 9x6x6 slab within the 15x6x6 remaining volume.
    # Orient blocks as 9x1x1.
    # Number of blocks = floor(9/9) * floor(6/1) * floor(6/1)
    blocks_3 = math.floor(9/9) * math.floor(6/1) * math.floor(6/1)

    # The final remaining volume is 6x6x6, where no 9x1x1 block can fit.

    # Total number of blocks is the sum of the blocks from each step.
    total_blocks = blocks_1 + blocks_2 + blocks_3

    # Print the final equation with the numbers calculated in each step.
    print(f"The largest number of blocks is the sum of blocks from each packing step:")
    print(f"Total blocks = {blocks_1} + {blocks_2} + {blocks_3} = {total_blocks}")

solve_packing()