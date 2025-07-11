import math

def solve_packing():
    """
    Calculates the maximum number of 9x1x1 blocks that fit inside a 15x15x15 cube.
    """
    # Stage 1: Pack a 9x15x15 section of the 15x15x15 cube.
    # We align the block's length of 9 with a cube edge of 15.
    # Number of blocks = floor(15/9) * floor(15/1) * floor(15/1)
    blocks_stage1 = (15 // 9) * (15 // 1) * (15 // 1)

    # This leaves a remaining space of (15-9) x 15 x 15 = 6x15x15.

    # Stage 2: Pack the 6x15x15 remaining space.
    # We must align the block's length of 9 with a 15-unit edge.
    # Number of blocks = floor(6/1) * floor(15/9) * floor(15/1)
    blocks_stage2 = (6 // 1) * (15 // 9) * (15 // 1)

    # This leaves a remaining space of 6 x (15-9) x 15 = 6x6x15.

    # Stage 3: Pack the 6x6x15 remaining space.
    # We must align the block's length of 9 with the 15-unit edge.
    # Number of blocks = floor(6/1) * floor(6/1) * floor(15/9)
    blocks_stage3 = (6 // 1) * (6 // 1) * (15 // 9)

    # The final remaining space is 6x6x(15-9) = 6x6x6, where no more blocks fit.

    total_blocks = blocks_stage1 + blocks_stage2 + blocks_stage3

    print("The maximum number of blocks is found by packing in stages:")
    print(f"1. In a 9x15x15 section: {blocks_stage1} blocks")
    print(f"2. In the remaining 6x15x15 space, we fill a 6x9x15 section with: {blocks_stage2} blocks")
    print(f"3. In the final 6x6x15 space, we fill a 6x6x9 section with: {blocks_stage3} blocks")
    print("\nThe total number of blocks is the sum of these stages:")
    print(f"{blocks_stage1} + {blocks_stage2} + {blocks_stage3} = {total_blocks}")

solve_packing()
<<<351>>>