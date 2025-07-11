import math

def solve():
    """
    Calculates the maximum number of 9x1x1 blocks that fit in a 15x15x15 cube.
    """
    cube_dim = 15
    block_len = 9
    block_wid = 1

    # Step 1: Pack a 15x15x15 space.
    # Orient the 9-unit length along one 15-unit edge.
    # The number of blocks is floor(15/9) * floor(15/1) * floor(15/1).
    num_blocks_1 = math.floor(cube_dim / block_len) * math.floor(cube_dim / block_wid) * math.floor(cube_dim / block_wid)
    
    # This leaves a remaining space of (15 - 9) x 15 x 15 = 6x15x15.
    rem_dim_1 = cube_dim - num_blocks_1 * block_len

    # Step 2: Pack the remaining 6x15x15 space.
    # We must orient the 9-unit length along a 15-unit edge.
    # The number of blocks is floor(6/1) * floor(15/9) * floor(15/1).
    num_blocks_2 = math.floor(rem_dim_1 / block_wid) * math.floor(cube_dim / block_len) * math.floor(cube_dim / block_wid)

    # This leaves a remaining space of 6 x (15 - 9) x 15 = 6x6x15.
    rem_dim_2 = cube_dim - math.floor(cube_dim / block_len) * block_len

    # Step 3: Pack the remaining 6x6x15 space.
    # We must orient the 9-unit length along the 15-unit edge.
    # The number of blocks is floor(6/1) * floor(6/1) * floor(15/9).
    num_blocks_3 = math.floor(rem_dim_1 / block_wid) * math.floor(rem_dim_2 / block_wid) * math.floor(cube_dim / block_len)

    # The final remaining space is 6x6x6, where no 9x1x1 block can fit.
    num_blocks_4 = 0

    # Step 4: Sum the blocks from all steps.
    total_blocks = num_blocks_1 + num_blocks_2 + num_blocks_3

    print("The calculation is done in stages:")
    print(f"1. In the 15x15x15 cube, we can fit {num_blocks_1} blocks. This leaves a 6x15x15 space.")
    print(f"2. In the remaining 6x15x15 space, we can fit {num_blocks_2} blocks. This leaves a 6x6x15 space.")
    print(f"3. In the remaining 6x6x15 space, we can fit {num_blocks_3} blocks. This leaves a 6x6x6 space.")
    print("4. No more blocks can fit in the final 6x6x6 space.")
    print("\nThe largest number of blocks is the sum of the blocks from each stage:")
    print(f"{num_blocks_1} + {num_blocks_2} + {num_blocks_3} = {total_blocks}")

solve()
<<<351>>>