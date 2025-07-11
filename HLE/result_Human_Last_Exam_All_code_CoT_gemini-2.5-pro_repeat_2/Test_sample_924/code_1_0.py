import math

def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes
    based on the derived mathematical and geometric constraints.
    """
    
    print("### Finding the Smallest Number of Green Cubes ###")
    
    # To find the minimum, we assume the core cube is Red (N_gcore = 0).
    # The formula for the total number of green cubes is:
    # N_total = N_gc + N_ge + N_gfc + N_gcore
    # From the face sum rule, 3*N_gc + 2*N_ge + N_gfc = 36, so N_gfc = 36 - 3*N_gc - 2*N_ge.
    # Substituting gives: N_total = 36 - 2*N_gc - N_ge + N_gcore.
    # With N_gcore = 0, we want to minimize (36 - 2*N_gc - N_ge),
    # which means we need to maximize the value of (2*N_gc + N_ge).

    max_val_to_maximize = -1
    best_config_min = {}

    # From geometric constraints, the number of green corners (N_gc) can only be 4, 5, 6, or 7.
    for n_gc in range(4, 8):
        # We need to find the maximum possible N_ge for this N_gc.
        # Constraints on N_ge arise from 0 <= N_gfc <= 6 and 0 <= N_ge <= 12.
        # 0 <= 36 - 3*n_gc - 2*n_ge <= 6
        # This simplifies to: (30 - 3*n_gc) / 2 <= n_ge <= (36 - 3*n_gc) / 2
        
        # We want the largest integer N_ge that satisfies the constraints.
        max_n_ge = math.floor((36 - 3 * n_gc) / 2)
        max_n_ge = min(max_n_ge, 12) # Can't have more than 12 edges.
        
        min_n_ge_check = math.ceil((30 - 3 * n_gc) / 2)
        
        if max_n_ge >= min_n_ge_check:
            current_val = 2 * n_gc + max_n_ge
            if current_val > max_val_to_maximize:
                max_val_to_maximize = current_val
                n_gfc = 36 - 3 * n_gc - 2 * max_n_ge
                best_config_min = {'N_gc': n_gc, 'N_ge': max_n_ge, 'N_gfc': n_gfc, 'N_gcore': 0}

    min_total_green = 36 - max_val_to_maximize
    
    print(f"The smallest possible number of green cubes is: {min_total_green}")
    min_cfg = best_config_min
    print(f"This can be achieved with a configuration of:")
    print(f"{min_cfg['N_gc']} green corners + {min_cfg['N_ge']} green edges + {min_cfg['N_gfc']} green face-centers + {min_cfg['N_gcore']} green core cubes.")
    print(f"Equation: {min_cfg['N_gc']} + {min_cfg['N_ge']} + {min_cfg['N_gfc']} + {min_cfg['N_gcore']} = {min_total_green}\n")


    print("### Finding the Largest Number of Green Cubes ###")
    
    # To find the maximum, we assume the core cube is Green (N_gcore = 1).
    # N_total = 37 - 2*N_gc - N_ge.
    # To maximize N_total, we must minimize the value of (2*N_gc + N_ge).
    
    min_val_to_minimize = float('inf')
    best_config_max = {}
    
    for n_gc in range(4, 8):
        # We need to find the minimum possible N_ge for this N_gc using the same constraints.
        min_n_ge = math.ceil((30 - 3 * n_gc) / 2)
        min_n_ge = max(min_n_ge, 0) # Can't have negative edges.

        max_n_ge_check = math.floor((36 - 3 * n_gc) / 2)

        if min_n_ge <= max_n_ge_check:
            current_val = 2 * n_gc + min_n_ge
            if current_val < min_val_to_minimize:
                min_val_to_minimize = current_val
                n_gfc = 36 - 3 * n_gc - 2 * min_n_ge
                best_config_max = {'N_gc': n_gc, 'N_ge': min_n_ge, 'N_gfc': n_gfc, 'N_gcore': 1}
    
    max_total_green = 37 - min_val_to_minimize

    print(f"The largest possible number of green cubes is: {max_total_green}")
    max_cfg = best_config_max
    print(f"This can be achieved with a configuration of:")
    print(f"{max_cfg['N_gc']} green corners + {max_cfg['N_ge']} green edges + {max_cfg['N_gfc']} green face-centers + {max_cfg['N_gcore']} green core cube.")
    print(f"Equation: {max_cfg['N_gc']} + {max_cfg['N_ge']} + {max_cfg['N_gfc']} + {max_cfg['N_gcore']} = {max_total_green}")

solve_cube_problem()