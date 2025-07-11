import math

def solve_packing_problem():
    """
    Calculates the maximum number of 9x1x1 blocks that fit in a 15x15x15 cube.
    """
    cube_edge = 15
    block_l, block_w, block_h = 9, 1, 1

    print("To find the maximum number of 9x1x1 blocks in a 15x15x15 cube, we pack them in stages.")
    print("-" * 30)

    # --- Stage 1 ---
    print("Step 1: Fill the cube with blocks in the first orientation.")
    # Align the 9-unit side with an edge of the cube
    n1_x = math.floor(cube_edge / block_l)
    n1_y = math.floor(cube_edge / block_w)
    n1_z = math.floor(cube_edge / block_h)
    blocks_1 = n1_x * n1_y * n1_z
    
    print(f"We align the 9-unit side of the blocks along one 15-unit edge.")
    print(f"Number of blocks fitting: floor({cube_edge}/{block_l}) * floor({cube_edge}/{block_w}) * floor({cube_edge}/{block_h}) = {n1_x} * {n1_y} * {n1_z} = {blocks_1}")
    
    # Calculate remaining space
    rem_dim1 = cube_edge - (n1_x * block_l)
    print(f"This occupies a {n1_x * block_l}x{cube_edge}x{cube_edge} volume, leaving a {rem_dim1}x{cube_edge}x{cube_edge} space empty.")
    print("-" * 30)

    # --- Stage 2 ---
    print("Step 2: Fill the remaining space with a different orientation.")
    # The remaining space is 6x15x15.
    space_dims_2 = [rem_dim1, cube_edge, cube_edge]
    print(f"Now we pack the {space_dims_2[0]}x{space_dims_2[1]}x{space_dims_2[2]} empty space.")
    print("We must orient the 9-unit side of the blocks along one of the 15-unit edges.")
    
    n2_x = math.floor(space_dims_2[0] / block_w)
    n2_y = math.floor(space_dims_2[1] / block_l)
    n2_z = math.floor(space_dims_2[2] / block_h)
    blocks_2 = n2_x * n2_y * n2_z
    print(f"Number of blocks fitting: floor({space_dims_2[0]}/{block_w}) * floor({space_dims_2[1]}/{block_l}) * floor({space_dims_2[2]}/{block_h}) = {n2_x} * {n2_y} * {n2_z} = {blocks_2}")

    # Calculate new remaining space
    rem_dim2 = space_dims_2[1] - (n2_y * block_l)
    print(f"This leaves a new empty space of {space_dims_2[0]}x{rem_dim2}x{space_dims_2[2]}.")
    print("-" * 30)

    # --- Stage 3 ---
    print("Step 3: Fill the new remaining space.")
    # The remaining space is 6x6x15
    space_dims_3 = [space_dims_2[0], rem_dim2, space_dims_2[2]]
    print(f"Now we pack the {space_dims_3[0]}x{space_dims_3[1]}x{space_dims_3[2]} empty space.")
    print("We orient the 9-unit side of the blocks along the final 15-unit edge.")
    
    n3_x = math.floor(space_dims_3[0] / block_w)
    n3_y = math.floor(space_dims_3[1] / block_h)
    n3_z = math.floor(space_dims_3[2] / block_l)
    blocks_3 = n3_x * n3_y * n3_z
    print(f"Number of blocks fitting: floor({space_dims_3[0]}/{block_w}) * floor({space_dims_3[1]}/{block_h}) * floor({space_dims_3[2]}/{block_l}) = {n3_x} * {n3_y} * {n3_z} = {blocks_3}")

    # Calculate final remaining space
    rem_dim3 = space_dims_3[2] - (n3_z * block_l)
    print(f"The final remaining empty space is {space_dims_3[0]}x{space_dims_3[1]}x{rem_dim3}.")
    print("No more 9x1x1 blocks can fit into this 6x6x6 space.")
    print("-" * 30)

    # --- Final Calculation ---
    total_blocks = blocks_1 + blocks_2 + blocks_3
    print("Final Calculation:")
    print("The total number of blocks is the sum from all steps:")
    print(f"{blocks_1} + {blocks_2} + {blocks_3} = {total_blocks}")

solve_packing_problem()