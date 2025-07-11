def get_obstruction_groups():
    """
    This function identifies the homotopy-theoretic groups that classify
    the obstructions for the two paths of bundle automorphisms to be homotopic.

    The variables 'n' and 'k' from the problem statement are used
    symbolically in the output strings, as no concrete values were given.
    """

    # Symbolic variables from the problem description
    n = "n"
    rank = "2k"

    # The obstruction groups are derived from a spectral sequence argument.
    # They are the constituent groups that form the fundamental group
    # of the space of bundle automorphisms.

    # 1. The relevant homology group of the base space X
    # The term from the spectral sequence is Hom(H_tilde_{n-1}(X), ...),
    # so we list the homology group itself. The (n-1) in H_{n-1}
    # is an equation involving the variable n.
    homology_group_X = f"H_tilde_{{{n}-1}}(X)"

    # 2. The first homotopy group of SO(2k)
    # This corresponds to the E_2^{0,1} term of the spectral sequence.
    # The 1 in pi_1 is a constant, and 2k is an equation.
    homotopy_group_1 = f"pi_1(SO({rank}))"

    # 3. The (n+1)-th homotopy group of SO(2k)
    # This corresponds to the E_2^{n, n+1} term.
    # The (n+1) is an equation involving the variable n.
    homotopy_group_n_plus_1 = f"pi_{{{n}+1}}(SO({rank}))"

    obstruction_groups = [
        homology_group_X,
        homotopy_group_1,
        homotopy_group_n_plus_1,
    ]

    print("The homotopy-theoretic obstructions are classified by a group constructed from the following list of groups:")
    # Print the final list of groups, as requested.
    print(obstruction_groups)

if __name__ == '__main__':
    get_obstruction_groups()