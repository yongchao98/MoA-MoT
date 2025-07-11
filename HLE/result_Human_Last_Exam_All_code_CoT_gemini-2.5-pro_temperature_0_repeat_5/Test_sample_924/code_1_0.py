def solve_cube_puzzle():
    """
    This function calculates and explains the smallest and largest possible
    number of green cubes based on the puzzle's rules.
    """

    # Total number of small cubes
    total_cubes = 27

    # There are 8 corners, 12 edges, 6 face-centers, and 1 core cube.

    # --- To find the SMALLEST number of green cubes ---
    # This corresponds to the MAXIMUM number of red cubes.
    # From our logical deduction, a valid configuration with the maximum
    # number of red cubes is:
    # 4 red corners, 0 red edges, 6 red face-centers, and 1 red core cube.
    r_c_max = 4
    r_e_max = 0
    r_f_max = 6
    r_i_max = 1 # Make the core red to maximize the total
    
    max_red_cubes = r_c_max + r_e_max + r_f_max + r_i_max
    
    # The number of green cubes is derived from the number of red cubes.
    g_c_min = 8 - r_c_max
    g_e_min = 12 - r_e_max
    g_f_min = 6 - r_f_max
    g_i_min = 1 - r_i_max
    
    min_green_cubes = g_c_min + g_e_min + g_f_min + g_i_min

    print("To find the smallest number of green cubes, we maximize the number of red cubes.")
    print(f"A valid configuration consists of {r_c_max} red corners, {r_e_max} red edges, {r_f_max} red face-centers, and {r_i_max} red core cube.")
    print(f"This means we have {g_c_min} green corners, {g_e_min} green edges, {g_f_min} green face-centers, and {g_i_min} green core cube.")
    print(f"Smallest number of green cubes = {g_c_min} + {g_e_min} + {g_f_min} + {g_i_min} = {min_green_cubes}")
    print("-" * 20)

    # --- To find the LARGEST number of green cubes ---
    # This corresponds to the MINIMUM number of red cubes.
    # From our logical deduction, a valid configuration with the minimum
    # number of red cubes is:
    # 6 red corners, 0 red edges, 0 red face-centers, and 0 red core cubes.
    r_c_min = 6
    r_e_min = 0
    r_f_min = 0
    r_i_min = 0 # Make the core green to minimize the total
    
    min_red_cubes = r_c_min + r_e_min + r_f_min + r_i_min

    # The number of green cubes is derived from the number of red cubes.
    g_c_max = 8 - r_c_min
    g_e_max = 12 - r_e_min
    g_f_max = 6 - r_f_min
    g_i_max = 1 - r_i_min

    max_green_cubes = g_c_max + g_e_max + g_f_max + g_i_max

    print("To find the largest number of green cubes, we minimize the number of red cubes.")
    print(f"A valid configuration consists of {r_c_min} red corners, {r_e_min} red edges, {r_f_min} red face-centers, and {r_i_min} red core cube.")
    print(f"This means we have {g_c_max} green corners, {g_e_max} green edges, {g_f_max} green face-centers, and {g_i_max} green core cube.")
    print(f"Largest number of green cubes = {g_c_max} + {g_e_max} + {g_f_max} + {g_i_max} = {max_green_cubes}")
    print("-" * 20)
    
    print(f"The smallest possible number of green cubes is {min_green_cubes}.")
    print(f"The largest possible number of green cubes is {max_green_cubes}.")


solve_cube_puzzle()