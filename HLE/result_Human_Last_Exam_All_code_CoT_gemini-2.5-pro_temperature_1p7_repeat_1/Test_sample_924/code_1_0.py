def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes
    based on the derived valid configurations.
    """
    
    # Calculation for the smallest number of green cubes
    # This corresponds to the maximum number of red cubes.
    # From our analysis, a valid configuration is r_c=4, r_e=0, r_f=6, r_i=1.
    r_c_max_R = 4
    r_e_max_R = 0
    r_f_max_R = 6
    r_i_max_R = 1
    max_red_cubes = r_c_max_R + r_e_max_R + r_f_max_R + r_i_max_R
    min_green_cubes = 27 - max_red_cubes

    # Calculation for the largest number of green cubes
    # A valid configuration found is g_c=8, g_e=3, g_f=6, g_i=1.
    g_c_max_G = 8
    g_e_max_G = 3
    g_f_max_G = 6
    g_i_max_G = 1
    max_green_cubes = g_c_max_G + g_e_max_G + g_f_max_G + g_i_max_G

    print("To find the smallest number of green cubes, we maximize the number of red cubes.")
    print(f"A valid configuration for red cubes is:")
    print(f"  - {r_c_max_R} red corners")
    print(f"  - {r_e_max_R} red edges")
    print(f"  - {r_f_max_R} red face-centers")
    print(f"  - {r_i_max_R} red core cube")
    print(f"Maximum number of red cubes = {r_c_max_R} + {r_e_max_R} + {r_f_max_R} + {r_i_max_R} = {max_red_cubes}")
    print(f"Smallest number of green cubes = 27 - {max_red_cubes} = {min_green_cubes}")
    print("-" * 30)
    print("To find the largest number of green cubes, we find a valid configuration for them directly.")
    print(f"A valid configuration for green cubes is:")
    print(f"  - {g_c_max_G} green corners")
    print(f"  - {g_e_max_G} green edges")
    print(f"  - {g_f_max_G} green face-centers")
    print(f"  - {g_i_max_G} green core cube")
    print(f"Largest number of green cubes = {g_c_max_G} + {g_e_max_G} + {g_f_max_G} + {g_i_max_G} = {max_green_cubes}")
    print("-" * 30)
    print(f"The smallest possible number of green cubes is: {min_green_cubes}")
    print(f"The largest possible number of green cubes is: {max_green_cubes}")

solve_cube_problem()