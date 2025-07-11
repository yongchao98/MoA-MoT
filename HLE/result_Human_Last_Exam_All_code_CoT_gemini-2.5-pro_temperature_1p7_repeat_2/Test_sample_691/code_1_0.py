def solve_topology_problem():
    """
    This function solves the topological problem by explaining the steps and printing the result.
    """

    # Step 1: Define the components.
    g_pants = 0  # Genus of a pair of pants (a sphere with 3 holes)
    n_pants = 3  # Boundary components of a pair of pants

    # Step 2: Calculate Euler characteristic of a single pair of pants.
    chi_pants = 2 - 2 * g_pants - n_pants
    print(f"Step 1: A pair of pants is a surface of genus g={g_pants} with n={n_pants} boundaries.")
    print(f"Its Euler characteristic chi = 2 - 2*g - n = 2 - 2*{g_pants} - {n_pants} = {chi_pants}.")
    print("-" * 20)

    # Step 3: Calculate Euler characteristic of the sewn-together space S'.
    # Gluing along two circles (chi=0) doesn't change the sum of characteristics.
    chi_S_prime = 2 * chi_pants
    print(f"Step 2: Sewing two pants along two leg openings gives a new surface S'.")
    print(f"The Euler characteristic of S' is the sum of the two pants' characteristics: {chi_pants} + {chi_pants} = {chi_S_prime}.")
    print("-" * 20)
    
    # Step 4: Determine the genus of S'.
    n_S_prime = 2  # Two remaining boundaries (the waistbands)
    # chi_S_prime = 2 - 2*g_S_prime - n_S_prime
    # -2 = 2 - 2*g - 2 => -2 = -2*g => g=1
    g_S_prime = (2 - n_S_prime - chi_S_prime) // 2
    print(f"Step 3: The surface S' has n={n_S_prime} boundaries (the two waistbands).")
    print(f"We find its genus g' using chi = 2 - 2*g' - n.")
    print(f"{chi_S_prime} = 2 - 2*g' - {n_S_prime}  =>  -2 = -2*g'  =>  g' = {g_S_prime}.")
    print("This means S' is a torus with two holes.")
    print("-" * 20)
    
    # Step 5: Alternative calculation for the final space X.
    print(f"Step 4: An alternative view simplifies finding the fundamental group.")
    print("  a) First, collapse the waistband of one pair of pants. A pair of pants is a disk with two holes.")
    print("     Collapsing its outer boundary gives a sphere with two punctures, which is a cylinder (S^1 x I).")
    print("  b) The problem becomes: take two cylinders and glue their corresponding circular ends together.")
    print("     (End 1 of Cyl 1 to End 1 of Cyl 2; End 2 of Cyl 1 to End 2 of Cyl 2).")
    print("  c) This construction results in a Torus (S^1 x S^1).")
    print("-" * 20)

    # Step 6: Final Answer.
    final_group_latex = "\\mathbb{Z} \\times \\mathbb{Z}"
    group_generator_1 = "\\mathbb{Z}"
    group_generator_2 = "\\mathbb{Z}"
    print(f"Step 5: The fundamental group of the torus (S^1 x S^1) is pi_1(S^1) x pi_1(S^1).")
    print(f"The final fundamental group is: {group_generator_1} x {group_generator_2}")
    
solve_topology_problem()