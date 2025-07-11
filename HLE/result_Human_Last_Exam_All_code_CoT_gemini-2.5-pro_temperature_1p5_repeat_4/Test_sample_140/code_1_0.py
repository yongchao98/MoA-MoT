import math

def solve_packing_problem():
    """
    Calculates the maximum number of 9x1x1 blocks that can fit
    inside a 15x15x15 cube using a step-by-step packing strategy.
    """
    
    # --- Stage 1: Initial packing in the 15x15x15 cube ---
    # We align the block's 9-unit side with one of the cube's axes.
    # This fills a 9x15x15 section of the cube.
    
    s1_x, s1_y, s1_z = 15, 15, 15
    b_l, b_w, b_h = 9, 1, 1
    
    # Calculate blocks in the first section
    # Fit along x-axis: math.floor(15/9) = 1
    # Fit along y-axis: math.floor(15/1) = 15
    # Fit along z-axis: math.floor(15/1) = 15
    blocks_stage1 = math.floor(s1_x / b_l) * math.floor(s1_y / b_w) * math.floor(s1_z / b_h)
    
    # The remaining space is (15-9) x 15 x 15 = 6x15x15
    
    # --- Stage 2: Packing the first remaining space (6x15x15) ---
    # We must orient the block's 9-unit side along a 15-unit edge.
    s2_x, s2_y, s2_z = 6, 15, 15
    
    # Calculate blocks in the second section
    # Fit along x-axis (6): math.floor(6/1) = 6
    # Fit along y-axis (15): math.floor(15/9) = 1
    # Fit along z-axis (15): math.floor(15/1) = 15
    blocks_stage2 = math.floor(s2_x / b_w) * math.floor(s2_y / b_l) * math.floor(s2_z / b_h)

    # The remaining space is 6 x (15-9) x 15 = 6x6x15
    
    # --- Stage 3: Packing the second remaining space (6x6x15) ---
    # We must orient the block's 9-unit side along the 15-unit edge.
    s3_x, s3_y, s3_z = 6, 6, 15

    # Calculate blocks in the third section
    # Fit along x-axis (6): math.floor(6/1) = 6
    # Fit along y-axis (6): math.floor(6/1) = 6
    # Fit along z-axis (15): math.floor(15/9) = 1
    blocks_stage3 = math.floor(s3_x / b_w) * math.floor(s3_y / b_h) * math.floor(s3_z / b_l)

    # The remaining space is 6x6x(15-9) = 6x6x6.
    # No 9x1x1 block can fit in a 6x6x6 cube.

    # --- Final Calculation ---
    total_blocks = blocks_stage1 + blocks_stage2 + blocks_stage3
    
    print("The problem is solved by packing the 15x15x15 cube in stages:")
    print(f"1. A 9x15x15 volume is filled with {blocks_stage1} blocks.")
    print(f"2. In the remaining 6x15x15 space, a 6x9x15 volume is filled with {blocks_stage2} blocks.")
    print(f"3. In the remaining 6x6x15 space, a 6x6x9 volume is filled with {blocks_stage3} blocks.")
    print("The final remaining space is a 6x6x6 cube, which cannot fit any more blocks.")
    print("\nThe final equation for the total number of blocks is:")
    print(f"{blocks_stage1} + {blocks_stage2} + {blocks_stage3} = {total_blocks}")

solve_packing_problem()