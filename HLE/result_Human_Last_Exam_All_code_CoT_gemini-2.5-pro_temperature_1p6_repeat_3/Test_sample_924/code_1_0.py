def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes
    in a 3x3x3 cube arrangement based on the given rules.
    """
    # Physical constraints on the number of cubes of each type
    max_Cc = 8  # Corner cubes
    max_Ce = 12 # Edge cubes
    max_Cf = 6  # Face-center cubes

    # Variables to store the results
    min_surface_green = float('inf')
    max_surface_green = float('-inf')
    min_config = {}
    max_config = {}

    # Iterate through all possible numbers of green surface cubes
    for Cc in range(max_Cc + 1):
        for Ce in range(max_Ce + 1):
            for Cf in range(max_Cf + 1):
                # Check if the combination satisfies the core equation for face patterns
                if 3 * Cc + 2 * Ce + 1 * Cf == 36:
                    surface_green_count = Cc + Ce + Cf
                    
                    # Check for minimum
                    if surface_green_count < min_surface_green:
                        min_surface_green = surface_green_count
                        min_config = {'Cc': Cc, 'Ce': Ce, 'Cf': Cf}
                        
                    # Check for maximum
                    if surface_green_count > max_surface_green:
                        max_surface_green = surface_green_count
                        max_config = {'Cc': Cc, 'Ce': Ce, 'Cf': Cf}

    # Smallest possible number: min surface green + a red core cube (0 green)
    min_total_green = min_surface_green + 0
    
    # Largest possible number: max surface green + a green core cube (1 green)
    max_total_green = max_surface_green + 1

    # --- Output the results ---
    
    # Smallest Number
    print(f"Smallest possible number of green cubes: {min_total_green}")
    min_cc = min_config['Cc']
    min_ce = min_config['Ce']
    min_cf = min_config['Cf']
    print(f"This is achieved with {min_cc} green corners, {min_ce} green edges, {min_cf} green face-centers, and 0 green core cubes.")
    print(f"Checking the face constraint equation: 3 * {min_cc} + 2 * {min_ce} + 1 * {min_cf} = {3*min_cc} + {2*min_ce} + {1*min_cf} = 36\n")

    # Largest Number
    print(f"Largest possible number of green cubes: {max_total_green}")
    max_cc = max_config['Cc']
    max_ce = max_config['Ce']
    max_cf = max_config['Cf']
    print(f"This is achieved with {max_cc} green corners, {max_ce} green edges, {max_cf} green face-centers, and 1 green core cube.")
    print(f"Checking the face constraint equation: 3 * {max_cc} + 2 * {max_ce} + 1 * {max_cf} = {3*max_cc} + {2*max_ce} + {1*max_cf} = 36")


solve_cube_problem()