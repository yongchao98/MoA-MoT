import math

def solve_packing_problem():
    """
    Calculates the maximum number of 9x1x1 blocks that can fit inside a 15x15x15 cube.
    """
    cube_dim = 15
    block_dims = [9, 1, 1]

    # --- Stage 1: Main Packing ---
    # Pack the 15x15x15 cube with blocks aligned as 9x1x1.
    # We can fit 1 block along the 9-unit dimension and 15 along the other two.
    # This fills a 9x15x15 volume.
    blocks_stage1 = (cube_dim // block_dims[0]) * (cube_dim // block_dims[1]) * (cube_dim // block_dims[2])
    
    # The remaining space is a 6x15x15 block.
    rem_space1 = [cube_dim % block_dims[0], cube_dim, cube_dim] # [6, 15, 15]

    # --- Stage 2: Packing the first remaining space ---
    # We must reorient the blocks to 1x9x1 to fit in the 6x15x15 space.
    blocks_stage2 = (rem_space1[0] // block_dims[1]) * (rem_space1[1] // block_dims[0]) * (rem_space1[2] // block_dims[2])
    
    # The remaining space from this stage is 6x6x15.
    rem_space2 = [rem_space1[0], rem_space1[1] % block_dims[0], rem_space1[2]] # [6, 6, 15]
    
    # --- Stage 3: Packing the second remaining space ---
    # We must reorient the blocks to 1x1x9 to fit in the 6x6x15 space.
    blocks_stage3 = (rem_space2[0] // block_dims[1]) * (rem_space2[1] // block_dims[2]) * (rem_space2[2] // block_dims[0])
    
    # The final remaining space is 6x6x6, where no 9x1x1 block can fit.

    # --- Final Calculation ---
    total_blocks = blocks_stage1 + blocks_stage2 + blocks_stage3
    
    print("This problem can be solved by packing in stages:")
    print(f"Stage 1: In the main 15x15x15 volume, we can fit {blocks_stage1} blocks.")
    print(f"This leaves a 6x15x15 space.")
    print(f"Stage 2: In the 6x15x15 space, we can fit another {blocks_stage2} blocks by turning them.")
    print(f"This leaves a 6x6x15 space.")
    print(f"Stage 3: In the 6x6x15 space, we can fit a final {blocks_stage3} blocks by turning them again.")
    print(f"The final remaining space is 6x6x6, which is too small for any more blocks.")
    print("\nTotal number of blocks:")
    print(f"{blocks_stage1} + {blocks_stage2} + {blocks_stage3} = {total_blocks}")

solve_packing_problem()