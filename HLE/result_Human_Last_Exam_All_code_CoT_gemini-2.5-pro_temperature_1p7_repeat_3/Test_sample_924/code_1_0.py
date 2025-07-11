def find_green_cubes_range():
    """
    This function finds the minimum and maximum possible number of green cubes
    by iterating through all possible combinations of green corner, edge, and
    face-center cubes that satisfy the problem's constraints.
    """
    
    # There are 8 corners, 12 edges, and 6 face-centers in a 3x3x3 cube.
    max_gc = 8
    max_ge = 12
    max_gf = 6

    min_total_green = float('inf')
    max_total_green = float('-inf')

    # We will store the configurations for min and max values.
    min_config = {}
    max_config = {}

    # Iterate through all possible numbers of green cubes of each type.
    # gc: green corners, ge: green edges, gf: green face-centers
    for gc in range(max_gc + 1):
        for ge in range(max_ge + 1):
            for gf in range(max_gf + 1):
                # The total number of green facelets must be 36.
                # A green corner contributes 3 green faces, an edge 2, and a face-center 1.
                if 3 * gc + 2 * ge + 1 * gf == 36:
                    
                    # Based on analysis, some configurations are geometrically impossible.
                    # Case 1: If all corners are green (gc=8), it forces all face-centers to be green,
                    # which is not a solution to the equation. So gc=8 is impossible.
                    if gc == 8:
                        continue
                    
                    # Case 2: If all face-centers are green (gf=6) and all edges are green (ge=12),
                    # the middle row/column on each face would have 3 green cubes, violating the rule.
                    if gf == 6 and ge == 12:
                        continue

                    # For the minimum number of green cubes, the single core cube should be red (g_core=0).
                    # N_g_min = gc + ge + gf + 0
                    current_sum = gc + ge + gf
                    if current_sum < min_total_green:
                        min_total_green = current_sum
                        # g_core is 0 for the minimum case
                        min_config = {'gc': gc, 'ge': ge, 'gf': gf, 'g_core': 0}

                    # For the maximum number of green cubes, the single core cube should be green (g_core=1).
                    # N_g_max = gc + ge + gf + 1
                    if current_sum > max_total_green:
                        max_total_green = current_sum
                        # g_core is 1 for the maximum case
                        max_config = {'gc': gc, 'ge': ge, 'gf': gf, 'g_core': 1}

    # The actual minimum requires g_core = 0
    final_min = min_total_green + 0
    
    # The actual maximum requires g_core = 1
    final_max = max_total_green + 1
    
    print("To find the smallest possible number of green cubes:")
    mc = min_config
    print(f"We can use {mc['gc']} green corners, {mc['ge']} green edges, {mc['gf']} green face-centers, and {mc['g_core']} green core cubes.")
    print("The number of green faces is:")
    print(f"3 * {mc['gc']} + 2 * {mc['ge']} + 1 * {mc['gf']} = {3*mc['gc'] + 2*mc['ge'] + 1*mc['gf']}")
    print("The total number of green cubes is:")
    print(f"{mc['gc']} + {mc['ge']} + {mc['gf']} + {mc['g_core']} = {final_min}")
    print("\nThe smallest possible number of green cubes is", final_min)
    
    print("\n" + "-"*50 + "\n")

    print("To find the largest possible number of green cubes:")
    xc = max_config
    print(f"We can use {xc['gc']} green corners, {xc['ge']} green edges, {xc['gf']} green face-centers, and {xc['g_core']} green core cubes.")
    print("The number of green faces is:")
    print(f"3 * {xc['gc']} + 2 * {xc['ge']} + 1 * {xc['gf']} = {3*xc['gc'] + 2*xc['ge'] + 1*xc['gf']}")
    print("The total number of green cubes is:")
    print(f"{xc['gc']} + {xc['ge']} + {xc['gf']} + {xc['g_core']} = {final_max}")
    print("\nThe largest possible number of green cubes is", final_max)


find_green_cubes_range()