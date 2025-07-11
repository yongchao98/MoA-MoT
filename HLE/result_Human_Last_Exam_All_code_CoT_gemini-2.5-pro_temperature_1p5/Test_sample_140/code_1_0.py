import math

def solve():
    """
    Calculates the largest number of 9x1x1 blocks that fit inside a 15x15x15 cube
    using a greedy slicing algorithm.
    """
    cube_dim = 15
    block_len = 9
    block_w = 1
    block_h = 1

    print(f"Starting with a cube of size {cube_dim}x{cube_dim}x{cube_dim} and blocks of size {block_len}x{block_w}x{block_h}.")
    print("-" * 30)

    # Step 1: Pack the first section
    print("Step 1: Pack a 9x15x15 section of the cube.")
    dim1 = block_len
    dim2 = cube_dim
    dim3 = cube_dim
    blocks_step1 = (dim1 // block_len) * (dim2 // block_w) * (dim3 // block_h)
    print(f"Number of blocks in this section: ({dim1}/{block_len}) * ({dim2}/{block_w}) * ({dim3}/{block_h}) = {blocks_step1}")
    remaining_dim1 = cube_dim - block_len
    print(f"Remaining space is a {remaining_dim1}x{cube_dim}x{cube_dim} cuboid.\n")


    # Step 2: Pack the remaining 6x15x15 space
    print("Step 2: Pack the remaining 6x15x15 section.")
    print("The 9-unit block length must be aligned with a 15-unit edge.")
    # Dimensions of this section are (remaining_dim1, cube_dim, block_len) = (6, 15, 9)
    dim1 = remaining_dim1
    dim2 = block_len
    dim3 = cube_dim
    blocks_step2 = (dim1 // block_w) * (dim2 // block_len) * (dim3 // block_h)
    print(f"Number of blocks in this section: ({dim1}/{block_w}) * ({dim2}/{block_len}) * ({dim3}/{block_h}) = {blocks_step2}")
    remaining_dim2 = cube_dim - block_len
    print(f"Remaining space is a {remaining_dim1}x{remaining_dim2}x{cube_dim} cuboid.\n")


    # Step 3: Pack the remaining 6x6x15 space
    print("Step 3: Pack the remaining 6x6x15 section.")
    print("The 9-unit block length must be aligned with the 15-unit edge.")
    # Dimensions of this section are (remaining_dim1, remaining_dim2, block_len) = (6, 6, 9)
    dim1 = remaining_dim1
    dim2 = remaining_dim2
    dim3 = block_len
    blocks_step3 = (dim1 // block_w) * (dim2 // block_h) * (dim3 // block_len)
    print(f"Number of blocks in this section: ({dim1}/{block_w}) * ({dim2}/{block_h}) * ({dim3}/{block_len}) = {blocks_step3}")
    remaining_dim3 = cube_dim - block_len
    print(f"Remaining space is a {remaining_dim1}x{remaining_dim2}x{remaining_dim3} cuboid.\n")

    # Step 4: Final remaining space
    print("Step 4: The final remaining space is 6x6x6.")
    blocks_step4 = 0
    print("A 9x1x1 block cannot fit into this space. No more blocks can be packed.\n")
    print("-" * 30)

    # Final Calculation
    total_blocks = blocks_step1 + blocks_step2 + blocks_step3
    print("The total number of blocks is the sum of the blocks from each step:")
    print(f"{blocks_step1} + {blocks_step2} + {blocks_step3} = {total_blocks}")

solve()
<<<351>>>