def print_obstruction_groups():
    """
    This function prints the names of the homology and homotopy groups that
    are the building blocks for the homotopy-theoretic obstructions.
    The variables 'n' and 'k' are used as placeholders for the dimensions
    given in the problem.
    """
    n_dim = "n"
    rank = "2k"

    homotopy_group = f"pi_{n_dim}+1(SO({rank}))"
    homology_group_1 = f"H_{n_dim}-1(X)"
    homology_group_2 = f"H_{n_dim}-2(X)"

    print("The homotopy-theoretic obstructions for the two paths to be homotopic are determined by the following groups:")
    print(f"1. A homotopy group of SO({rank}): {homotopy_group}")
    print(f"2. A homology group of X: {homology_group_1}")
    print(f"3. Another homology group of X: {homology_group_2}")

print_obstruction_groups()