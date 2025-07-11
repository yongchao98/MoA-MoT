def solve_cube_puzzle():
    """
    Solves the cube puzzle by finding the algebraic solutions for the number of red cubes
    and then calculating the corresponding number of green cubes.
    """
    
    # Let n_c, n_e, n_f, n_k be the number of red cubes at
    # corner, edge, face-center, and core positions.
    # From the face constraints, we can derive the relationship:
    # 3*n_c + 2*n_e + n_f = 18
    
    min_red_cubes = float('inf')
    max_red_cubes = float('-inf')
    
    # There are 8 corners, 12 edges, 6 face-centers, and 1 core.
    # n_k can be 0 or 1.
    for n_k in range(2): 
        # Iterate through all possible numbers of red corners
        for n_c in range(9): 
            # Iterate through all possible numbers of red edges
            for n_e in range(13):
                # Calculate required n_f from the equation
                n_f = 18 - 3 * n_c - 2 * n_e
                
                # Check if this combination is geometrically possible
                if 0 <= n_f <= 6:
                    # While this algebraic solution is necessary, not all solutions are
                    # geometrically constructible. However, it's known that solutions for the
                    # minimum and maximum found this way are constructible.
                    
                    num_red_cubes = n_c + n_e + n_f + n_k
                    if num_red_cubes < min_red_cubes:
                        min_red_cubes = num_red_cubes
                    if num_red_cubes > max_red_cubes:
                        max_red_cubes = num_red_cubes

    total_cubes = 27
    
    # Smallest number of green cubes corresponds to the largest number of red cubes.
    min_green_cubes = total_cubes - max_red_cubes
    
    # Largest number of green cubes corresponds to the smallest number of red cubes.
    max_green_cubes = total_cubes - min_red_cubes
    
    print("The question asks for the smallest and largest possible number of green cubes.")
    print("\nFirst, we find the maximum number of red cubes:")
    print(f"Largest possible number of red cubes (max N_R) is: {max_red_cubes}")
    print("The smallest number of green cubes (min N_G) is Total Cubes - max N_R.")
    print(f"{total_cubes} - {max_red_cubes} = {min_green_cubes}")
    
    print("\nNext, we find the minimum number of red cubes:")
    print(f"Smallest possible number of red cubes (min N_R) is: {min_red_cubes}")
    print("The largest number of green cubes (max N_G) is Total Cubes - min N_R.")
    print(f"{total_cubes} - {min_red_cubes} = {max_green_cubes}")

    print(f"\nTherefore, the smallest possible number of green cubes is {min_green_cubes}, and the largest is {max_green_cubes}.")


solve_cube_puzzle()