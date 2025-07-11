def solve_genus_problem():
    """
    Solves the topology problem by calculating the smallest possible genus
    for a bounding surface.
    """
    # Given information about the initial surface Sigma
    g_sigma = 10  # Genus of Sigma
    b_sigma = 1   # Number of boundary components of Sigma

    # --- Part 1: Finding an upper bound by construction ---
    print("Part 1: Constructing a bounding surface Sigma' and finding its genus.")
    print("-" * 60)

    # Step 1: Calculate the Euler characteristic of Sigma.
    # The formula for a surface with genus g and b boundary components is chi = 2 - 2g - b.
    chi_sigma = 2 - 2 * g_sigma - b_sigma
    print("Step 1: Calculate the Euler characteristic of the initial surface Sigma.")
    print(f"Sigma has genus g = {g_sigma} and b = {b_sigma} boundary component.")
    print(f"The formula is: chi(Sigma) = 2 - 2*g - b")
    print(f"chi(Sigma) = 2 - 2*{g_sigma} - {b_sigma} = {chi_sigma}")
    print()

    # Step 2: Form a closed surface Sigma' by capping the boundary with a disk.
    # A disk has genus g=0, b=1, and its Euler characteristic is chi(Disk) = 1.
    chi_disk = 1
    print("Step 2: Form a closed surface Sigma' by capping the boundary.")
    print("Since the boundary is an unknotted circle, we can glue a disk D to it.")
    print(f"The Euler characteristic of a disk is chi(D) = {chi_disk}.")
    print()
    
    # Step 3: Calculate the Euler characteristic of the new closed surface Sigma'.
    # The Euler characteristic of the union is the sum of the individual characteristics.
    chi_sigma_prime = chi_sigma + chi_disk
    print("Step 3: Calculate the Euler characteristic of the new closed surface Sigma'.")
    print("The formula is: chi(Sigma') = chi(Sigma) + chi(D)")
    print(f"chi(Sigma') = {chi_sigma} + {chi_disk} = {chi_sigma_prime}")
    print()

    # Step 4: Calculate the genus of the new closed surface Sigma' from its Euler characteristic.
    # For a closed surface (b=0), the formula is chi = 2 - 2g.
    # Rearranging gives: 2g = 2 - chi  =>  g = (2 - chi) / 2.
    numerator = 2 - chi_sigma_prime
    g_sigma_prime = numerator // 2
    
    print("Step 4: Calculate the genus g' of the closed surface Sigma'.")
    print("For a closed surface, the formula is chi = 2 - 2g'.")
    print("Rearranging for g': 2g' = 2 - chi'")
    print(f"2 * g' = 2 - ({chi_sigma_prime})")
    print(f"2 * g' = {numerator}")
    print(f"g' = {numerator} / 2 = {g_sigma_prime}")
    print()
    print(f"This construction shows that for any such Sigma, there exists a containing surface Sigma' of genus {g_sigma_prime}.")
    print("-" * 60)
    print()

    # --- Part 2: Finding the lower bound ---
    print("Part 2: Proving this is the smallest possible genus.")
    print("-" * 60)
    print("The genus of a subsurface cannot be greater than the genus of the surface containing it.")
    print("Since Sigma is a subsurface of Sigma' (Sigma is contained in Sigma'), we must have:")
    print("g(Sigma') >= g(Sigma)")
    print(f"g(Sigma') >= {g_sigma}")
    print()
    print(f"This shows that any containing surface must have a genus of at least {g_sigma}.")
    print("-" * 60)
    print()

    # --- Conclusion ---
    print("Conclusion:")
    print(f"From Part 1, we know it is always possible to find a containing surface of genus {g_sigma_prime}.")
    print(f"From Part 2, we know the genus of any such surface cannot be less than {g_sigma}.")
    print(f"Since these bounds meet, the smallest possible positive integer for the genus g is {g_sigma_prime}.")
    
    return g_sigma_prime

# Execute the function to print the solution steps.
final_answer = solve_genus_problem()
print(f"<<<{final_answer}>>>")