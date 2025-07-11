def solve_cube_puzzle():
    """
    Calculates the minimum and maximum number of green cubes required to build a 3x3x3 cube
    under the rule that each row and column on every face has two green cubes and one red cube.
    """
    min_green_cubes = float('inf')
    max_green_cubes = float('-inf')
    
    min_config = {}
    max_config = {}

    # Total number of cubelets of each type
    total_corners = 8
    total_edges = 12
    total_face_centers = 6
    
    # Derived constraint for the number of green corners: 4 <= g_c <= 6
    for g_c in range(4, 6 + 1):
        # Iterate through possible numbers of green face-centers
        for g_f in range(0, total_face_centers + 1):
            
            # From the main equation: 3*g_c + 2*g_e + 1*g_f = 36
            # We can find the required number of green edges, g_e.
            # 2*g_e = 36 - 3*g_c - g_f
            
            numerator = 36 - 3 * g_c - g_f
            
            # g_e must be an integer, so the numerator must be even.
            if numerator % 2 == 0:
                g_e = numerator // 2
                
                # g_e must be a valid number of cubes.
                if 0 <= g_e <= total_edges:
                    
                    # This combination (g_c, g_e, g_f) is valid.
                    # Now, calculate the total green cubes, G = g_c + g_e + g_f + g_core.
                    # The core cube can be red (g_core=0) or green (g_core=1).
                    
                    # Case 1: Core cube is red (g_core = 0)
                    g_core_0 = 0
                    total_green_0 = g_c + g_e + g_f + g_core_0
                    
                    if total_green_0 < min_green_cubes:
                        min_green_cubes = total_green_0
                        min_config = {'g_c': g_c, 'g_e': g_e, 'g_f': g_f, 'g_core': g_core_0}
                    
                    if total_green_0 > max_green_cubes:
                        max_green_cubes = total_green_0
                        max_config = {'g_c': g_c, 'g_e': g_e, 'g_f': g_f, 'g_core': g_core_0}

                    # Case 2: Core cube is green (g_core = 1)
                    g_core_1 = 1
                    total_green_1 = g_c + g_e + g_f + g_core_1

                    if total_green_1 < min_green_cubes:
                        min_green_cubes = total_green_1
                        min_config = {'g_c': g_c, 'g_e': g_e, 'g_f': g_f, 'g_core': g_core_1}

                    if total_green_1 > max_green_cubes:
                        max_green_cubes = total_green_1
                        max_config = {'g_c': g_c, 'g_e': g_e, 'g_f': g_f, 'g_core': g_core_1}

    print("--- Smallest Possible Number of Green Cubes ---")
    min_c = min_config['g_c']
    min_e = min_config['g_e']
    min_f = min_config['g_f']
    min_core = min_config['g_core']
    print(f"To achieve the minimum of {min_green_cubes} green cubes:")
    print(f"Number of green corner cubes (g_c): {min_c}")
    print(f"Number of green edge cubes (g_e): {min_e}")
    print(f"Number of green face-center cubes (g_f): {min_f}")
    print(f"Number of green core cubes (g_core): {min_core}")
    print(f"Equation: {min_c} + {min_e} + {min_f} + {min_core} = {min_green_cubes}")
    
    print("\n--- Largest Possible Number of Green Cubes ---")
    max_c = max_config['g_c']
    max_e = max_config['g_e']
    max_f = max_config['g_f']
    max_core = max_config['g_core']
    print(f"To achieve the maximum of {max_green_cubes} green cubes:")
    print(f"Number of green corner cubes (g_c): {max_c}")
    print(f"Number of green edge cubes (g_e): {max_e}")
    print(f"Number of green face-center cubes (g_f): {max_f}")
    print(f"Number of green core cubes (g_core): {max_core}")
    print(f"Equation: {max_c} + {max_e} + {max_f} + {max_core} = {max_green_cubes}")

solve_cube_puzzle()