def solve_packing_problem():
    """
    Calculates the largest number of 9x1x1 blocks that fit inside a 15x15x15 cube
    by breaking the problem down into packing stages.
    """
    
    # Stage 1: Pack a 9x15x15 volume from the initial 15x15x15 cube.
    # The number of blocks is (15/1) * (15/1) * (9/9).
    term1 = 15 * 15

    # Stage 2: Pack the remaining 6x15x15 volume.
    # We pack a 6x15x9 section. The number of blocks is (6/1) * (15/1) * (9/9).
    term2 = 6 * 15

    # Stage 3: Pack the remaining 6x6x15 volume.
    # We pack a 6x6x9 section. The number of blocks is (6/1) * (6/1) * (9/9).
    term3 = 6 * 6
    
    # The final remaining 6x6x6 volume cannot be packed.
    
    total_blocks = term1 + term2 + term3

    print("The maximum number of blocks is found by summing the blocks from each packing stage.")
    print("The equation for the total is:")
    
    # Print the final equation with each number and the result.
    print(f"{term1} + {term2} + {term3} = {total_blocks}")

solve_packing_problem()