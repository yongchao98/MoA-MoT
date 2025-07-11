def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes
    based on the given face pattern rules.
    """
    
    possible_surface_greens = []
    
    # N_A is the number of faces with Type A pattern (2 green corners, 4 green edges, 0 green centers)
    # The number of faces with Type B pattern (3g_c, 2g_e, 1g_f) is N_B = 6 - N_A
    
    # Iterate through possible values for N_A. For C_g to be an integer, (18-N_A) must be divisible by 3.
    # This means N_A must be a multiple of 3.
    for N_A in range(0, 7, 3):
        N_B = 6 - N_A
        
        # Calculate C_g (number of green corner cubes)
        # 3 * C_g = 2 * N_A + 3 * N_B = 2 * N_A + 3 * (6 - N_A) = 18 - N_A
        C_g = (18 - N_A) / 3
        
        # Calculate E_g (number of green edge cubes)
        # 2 * E_g = 4 * N_A + 2 * N_B = 4 * N_A + 2 * (6 - N_A) = 2 * N_A + 12
        E_g = (2 * N_A + 12) / 2
        
        # Calculate F_g (number of green face-center cubes)
        # F_g = N_B = 6 - N_A
        F_g = 6 - N_A
        
        # All calculated cube counts must be integers. The loops ensure this.
        # Check if the cube counts are within the possible physical limits.
        # 8 corners, 12 edges, 6 face-centers in total.
        if 0 <= C_g <= 8 and 0 <= E_g <= 12 and 0 <= F_g <= 6:
            G_surface = C_g + E_g + F_g
            possible_surface_greens.append({
                "G_surface": G_surface,
                "C_g": int(C_g),
                "E_g": int(E_g),
                "F_g": int(F_g),
            })
            
    # Find the min and max surface green counts
    min_surface_config = min(possible_surface_greens, key=lambda x: x["G_surface"])
    max_surface_config = max(possible_surface_greens, key=lambda x: x["G_surface"])
    
    # Min green cubes = min surface greens + 0 (red core cube)
    min_total_greens = min_surface_config["G_surface"]
    min_C_g = min_surface_config["C_g"]
    min_E_g = min_surface_config["E_g"]
    min_F_g = min_surface_config["F_g"]
    min_I_g = 0 # Core cube is red for the minimum case
    
    # Max green cubes = max surface greens + 1 (green core cube)
    max_total_greens = max_surface_config["G_surface"] + 1
    max_C_g = max_surface_config["C_g"]
    max_E_g = max_surface_config["E_g"]
    max_F_g = max_surface_config["F_g"]
    max_I_g = 1 # Core cube is green for the maximum case
    
    print(f"The smallest possible number of green cubes is {int(min_total_greens)}.")
    print(f"This is achieved with {min_C_g} green corners, {min_E_g} green edges, {min_F_g} green face-centers, and {min_I_g} green core cube.")
    print(f"Equation: {int(min_total_greens)} = {min_C_g} + {min_E_g} + {min_F_g} + {min_I_g}")
    print("")
    print(f"The largest possible number of green cubes is {int(max_total_greens)}.")
    print(f"This is achieved with {max_C_g} green corners, {max_E_g} green edges, {max_F_g} green face-centers, and {max_I_g} green core cube.")
    print(f"Equation: {int(max_total_greens)} = {max_C_g} + {max_E_g} + {max_F_g} + {max_I_g}")

solve_cube_problem()