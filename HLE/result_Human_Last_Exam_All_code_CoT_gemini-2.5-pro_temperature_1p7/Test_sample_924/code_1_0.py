def solve_cube_puzzle():
    """
    Calculates the smallest and largest possible number of green cubes based on the derived formulas.
    The total number of red cubes (N_R) is given by:
    N_R = N_R_core + N_R_c + 6
    where N_R_core is 0 or 1, and the number of red corners N_R_c can be 2, 3, or 4.
    The number of green cubes N_G = 27 - N_R.
    """
    
    # Range for the number of red corners
    min_red_corners = 2
    max_red_corners = 4
    
    # To find the maximum number of green cubes, we need the minimum number of red cubes.
    # This happens when the core is green (N_R_core = 0) and red corners are minimized.
    min_red_core = 0
    min_N_R = min_red_core + min_red_corners + 6
    max_green_cubes = 27 - min_N_R
    
    # To find the smallest number of green cubes, we need the maximum number of red cubes.
    # This happens when the core is red (N_R_core = 1) and red corners are maximized.
    max_red_core = 1
    max_N_R = max_red_core + max_red_corners + 6
    min_green_cubes = 27 - max_N_R
    
    print("To find the smallest number of green cubes:")
    print(f"We maximize the number of red cubes. This requires the maximum number of red corners ({max_red_corners}) and a red core cube ({max_red_core}).")
    print(f"Maximum Red Cubes = {max_red_core} (core) + {max_red_corners} (corners) + 6 = {max_N_R}")
    print(f"Smallest Green Cubes = 27 - {max_N_R} = {min_green_cubes}")
    print("-" * 30)
    print("To find the largest number of green cubes:")
    print(f"We minimize the number of red cubes. This requires the minimum number of red corners ({min_red_corners}) and a green core cube ({min_red_core}).")
    print(f"Minimum Red Cubes = {min_red_core} (core) + {min_red_corners} (corners) + 6 = {min_N_R}")
    print(f"Largest Green Cubes = 27 - {min_N_R} = {max_green_cubes}")
    print("-" * 30)
    print(f"The smallest possible number of green cubes is {min_green_cubes}.")
    print(f"The largest possible number of green cubes is {max_green_cubes}.")

solve_cube_puzzle()