def print_obstruction_groups():
    """
    This function prints the list of groups that form the building blocks
    for the homotopy-theoretic obstructions.
    The variables n, k, and X are symbolic, as given in the problem statement.
    """
    n = "n"
    k = "k"
    X = "X"
    q = "q"

    print("The homotopy-theoretic obstructions for the paths phi_t and psi_t are classified by an abelian group that is an extension of two smaller groups.")
    print("The building blocks for this obstruction group are given by the following list of groups:")
    print("")

    # 1. The local obstruction component
    group1 = f"pi_1(SO(2*{k}))"
    print(f"1. A subgroup of {group1}, related to the 'local' difference between the paths at a point.")
    
    # 2. Part of the global obstruction component
    group2 = f"pi_2(SO(2*{k}))"
    print(f"2. The group {group2}, which affects the structure of the 'global' obstruction component.")
    
    # 3. The main global obstruction component related to the base space topology
    print(f"3. A series of groups coming from the topology of the base space {X}:")
    homology_group = f"H^({n}-1)({X}; pi_{q}(SO(2*{k})))"
    print(f"   The integral homology groups of X with coefficients in the homotopy groups of SO(2k), of the form {homology_group} for {q} = 1, 2, 3, ...")


print_obstruction_groups()