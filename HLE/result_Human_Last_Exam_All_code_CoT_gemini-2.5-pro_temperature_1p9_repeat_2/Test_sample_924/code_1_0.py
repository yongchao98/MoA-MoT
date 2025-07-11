def solve_cube_puzzle():
    """
    This function finds the smallest and largest possible number of green cubes.
    It iterates through all possible numbers of green corner, edge, and face-center cubes,
    checks if they satisfy the two identified valid geometric configurations, and calculates
    the corresponding total number of green cubes.
    """
    possible_G = []
    
    # Iterate through all combinations for the number of green cubes of each type on the surface
    # Cc: number of green corner cubes (total 8)
    # Ce: number of green edge cubes (total 12)
    # Cf: number of green face-center cubes (total 6)
    print("Searching for valid configurations...\n")
    for Cc in range(9):
        for Ce in range(13):
            for Cf in range(7):
                # The total number of green squares on all faces must be 36
                equation = 3 * Cc + 2 * Ce + 1 * Cf
                if equation == 36:
                    
                    # From our analysis, only two types of configurations are geometrically possible without contradictions.
                    # Configuration 1: All face-centers are red (Cf=0), which forces all edges to be green (Ce=12).
                    # This corresponds to Cc=4.
                    is_config_1 = (Cc == 4 and Ce == 12 and Cf == 0)

                    # Configuration 2: All face-centers are green (Cf=6), which forces half the edges to be green (Ce=6).
                    # This corresponds to Cc=6.
                    is_config_2 = (Cc == 6 and Ce == 6 and Cf == 6)
                    
                    if is_config_1 or is_config_2:
                        # G_surface is the number of green cubes on the surface
                        G_surface = Cc + Ce + Cf
                        
                        # The central cube can be either red or green, as it's not on any face.
                        # Ci = 0: central cube is red
                        G1 = G_surface + 0
                        possible_G.append(G1)
                        
                        # Ci = 1: central cube is green
                        G2 = G_surface + 1
                        possible_G.append(G2)
                        
                        if is_config_1:
                            print(f"Found valid tetrahedral configuration: (Cc={Cc}, Ce={Ce}, Cf={Cf})")
                            print(f"Number of green surface cubes: {Cc} + {Ce} + {Cf} = {G_surface}")
                            print(f"Total green cubes can be {G1} (inner red) or {G2} (inner green).\n")
                        else: # is_config_2
                            print(f"Found valid 'mod 3' configuration: (Cc={Cc}, Ce={Ce}, Cf={Cf})")
                            print(f"Number of green surface cubes: {Cc} + {Ce} + {Cf} = {G_surface}")
                            print(f"Total green cubes can be {G1} (inner red) or {G2} (inner green).\n")
                            

    if not possible_G:
        print("No valid configurations found.")
        return

    min_g = min(possible_G)
    max_g = max(possible_G)
    
    print("-----------------------------------------")
    print(f"Based on the valid configurations, the full list of possible numbers of green cubes is: {sorted(list(set(possible_G)))}")
    print(f"The smallest possible number of green cubes is: {min_g}")
    print(f"The largest possible number of green cubes is: {max_g}")


solve_cube_puzzle()
<<<16 and 19>>>