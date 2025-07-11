def calculate_symmetry_breaking(initial_group_n, residual_su_n):
    """
    Calculates and prints the details of a symmetry breaking from
    SU(N) to SU(M) x U(1).
    """

    def su_n_generators(n):
        """Calculates the number of generators for the SU(N) group."""
        return n**2 - 1

    def u_1_generators():
        """Returns the number of generators for the U(1) group."""
        return 1

    # Step 1: Analyze the initial group G = SU(N)
    G_name = f"SU({initial_group_n})"
    num_total_generators = su_n_generators(initial_group_n)
    print(f"The initial symmetry group is {G_name}.")
    print(f"The number of total generators for {G_name} is {initial_group_n}^2 - 1 = {num_total_generators}.")
    print("This confirms statement D is correct.\n")

    # Step 2: Analyze the residual group H = SU(M) x U(1)
    H_name = f"SU({residual_su_n}) x U(1)"
    num_su_residual_generators = su_n_generators(residual_su_n)
    num_u1_residual_generators = u_1_generators()
    num_unbroken_generators = num_su_residual_generators + num_u1_residual_generators
    print(f"The residual (unbroken) symmetry group is {H_name}.")
    print(f"The SU({residual_su_n}) factor has {residual_su_n}^2 - 1 = {num_su_residual_generators} generators (Confirms H).")
    print(f"The U(1) factor has {num_u1_residual_generators} generator (Confirms G).")
    print(f"The total number of unbroken generators is {num_su_residual_generators} + {num_u1_residual_generators} = {num_unbroken_generators} (Confirms I).\n")

    # Step 3: Calculate the number of broken generators
    num_broken_generators = num_total_generators - num_unbroken_generators
    print("The vacuum degeneracy is characterized by the number of broken generators.")
    print("This is the dimension of the initial group minus the dimension of the residual group.")
    print(f"Number of broken generators = {num_total_generators} - {num_unbroken_generators} = {num_broken_generators}\n")
    print("Statement B (Five broken generators) is incorrect.\n")

    # Step 4: Relate to physical particles via the Higgs mechanism
    # The number of massive gauge bosons equals the number of broken generators.
    num_massive_gauge_bosons = num_broken_generators
    print("In a non-Abelian gauge theory, the Higgs mechanism gives mass to gauge bosons.")
    print("The number of massive gauge bosons is equal to the number of broken generators.")
    print(f"Number of massive gauge bosons = {num_broken_generators}\n")
    print("This shows statement E is correct.\n")

    print("Conclusion:")
    print("The 'vacuum degeneracy condition' refers to the consequences of the vacuum state having a lower symmetry than the underlying theory.")
    print("The dimension of the space of degenerate vacua is equal to the number of broken generators (4 in this case).")
    print("In a gauge theory, this directly results in an equal number of massive gauge bosons.")
    print("While several options are true statements about the groups, option E describes the key physical outcome that quantifies this degeneracy.")

# Run the calculation for SU(3) -> SU(2) x U(1)
calculate_symmetry_breaking(initial_group_n=3, residual_su_n=2)