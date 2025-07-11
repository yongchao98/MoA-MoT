import math

def calculate_blocks_in_cuboid(cuboid_dims, block_dims):
    """
    Calculates the maximum number of blocks that can be packed into a cuboid.
    It checks all 6 possible orientations of the block.
    """
    from itertools import permutations
    
    max_blocks = 0
    # Check every permutation of the block's dimensions
    for p in set(permutations(block_dims)):
        if p[0] <= cuboid_dims[0] and p[1] <= cuboid_dims[1] and p[2] <= cuboid_dims[2]:
            num_blocks = (math.floor(cuboid_dims[0] / p[0]) *
                          math.floor(cuboid_dims[1] / p[1]) *
                          math.floor(cuboid_dims[2] / p[2]))
            if num_blocks > max_blocks:
                max_blocks = num_blocks
    return max_blocks

def solve_packing_problem():
    """
    Solves the problem of packing 9x1x1 blocks into a 15x15x15 cube.
    """
    cube_dim = 15
    block_dims = (9, 1, 1)

    print(f"Finding the largest number of {block_dims[0]}x{block_dims[1]}x{block_dims[2]} blocks in a {cube_dim}x{cube_dim}x{cube_dim} cube.")
    print("\nStep 1: Decompose the cube based on the block's largest dimension.")
    large_part = block_dims[0]
    small_part = cube_dim % large_part
    print(f"The cube's edge length {cube_dim} is decomposed into: {cube_dim} = 1 * {large_part} + {small_part}")

    # These are the sub-cuboids resulting from the (9+6)^3 decomposition, and their counts
    sub_cuboids = {
        (9, 9, 9): 1,
        (9, 9, 6): 3,
        (9, 6, 6): 3,
        (6, 6, 6): 1
    }

    print("\nStep 2: Calculate how many blocks fit into each type of sub-cuboid.")
    total_blocks = 0
    calculation_parts = []
    
    for dims, count in sub_cuboids.items():
        blocks_in_one = calculate_blocks_in_cuboid(dims, block_dims)
        total_for_type = blocks_in_one * count
        total_blocks += total_for_type
        print(f"- For the {count} cuboid(s) of size {dims[0]}x{dims[1]}x{dims[2]}:")
        print(f"  {blocks_in_one} blocks fit in each.")
        calculation_parts.append(f"{blocks_in_one} * {count}")

    print("\nStep 3: Sum the blocks from all sub-cuboids.")
    
    final_equation = " + ".join(calculation_parts)
    
    # To show the intermediate products before the final sum
    intermediate_sums = []
    for part in calculation_parts:
        # Evaluate "num * count"
        intermediate_sums.append(str(eval(part)))
        
    print(f"Total blocks = {final_equation}")
    print(f"             = {' + '.join(intermediate_sums)}")
    print(f"             = {total_blocks}")
    
    print(f"\nThe largest number of blocks is {total_blocks}.")

solve_packing_problem()
<<<351>>>