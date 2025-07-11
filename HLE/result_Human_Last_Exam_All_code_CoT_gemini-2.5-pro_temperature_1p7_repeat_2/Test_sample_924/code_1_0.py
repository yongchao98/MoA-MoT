def solve_cube_puzzle():
    """
    Calculates the smallest and largest possible number of green cubes
    based on the derived formula G = g_c + 12 + g_m.
    """
    print("The total number of green cubes (G) is determined by the number of green corners (g_c) and green central cubes (g_m).")
    print("The derived formula is: G = g_c + 12 + g_m")
    print("-" * 20)

    # Calculate the minimum number of green cubes
    print("To find the smallest number of green cubes, we use the minimum possible values:")
    g_c_min = 4
    g_m_min = 0
    print(f"Minimum green corners (g_c): {g_c_min}")
    print(f"Minimum green central cubes (g_m): {g_m_min}")
    
    g_min = g_c_min + 12 + g_m_min
    
    print(f"Smallest G = {g_c_min} + 12 + {g_m_min}")
    print(f"The smallest possible number of green cubes is: {g_min}")
    print("-" * 20)

    # Calculate the maximum number of green cubes
    print("To find the largest number of green cubes, we use the maximum possible values:")
    g_c_max = 6
    g_m_max = 1
    print(f"Maximum green corners (g_c): {g_c_max}")
    print(f"Maximum green central cubes (g_m): {g_m_max}")

    g_max = g_c_max + 12 + g_m_max
    
    print(f"Largest G = {g_c_max} + 12 + {g_m_max}")
    print(f"The largest possible number of green cubes is: {g_max}")
    print("-" * 20)

solve_cube_puzzle()