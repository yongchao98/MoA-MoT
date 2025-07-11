def solve_group_problem():
    """
    This function provides the solution to the question about the number of
    finite groups with a maximal product-free set of size 2.

    The solution is based on a known classification theorem from finite group theory.
    """

    # According to the classification theorem by V. A. Belonogov,
    # a finite group G has a maximal product-free set of size 2
    # if and only if it is isomorphic to one of the following groups.
    groups_list = [
        "Z_4 (The cyclic group of order 4)",
        "Z_2 x Z_2 (The Klein four-group, order 4)",
        "Z_5 (The cyclic group of order 5)",
        "Z_6 (The cyclic group of order 6)",
        "D_10 (The dihedral group of order 10)",
        "A_4 (The alternating group on 4 elements, order 12)",
        "D_12 (The dihedral group of order 12)"
    ]

    # These 7 groups are non-isomorphic, so they are all distinct for the purpose
    # of this problem.

    number_of_groups = len(groups_list)

    print("The finite groups that contain a maximal by inclusion product-free set of size 2 are:")
    for group in groups_list:
        print(f"- {group}")

    print("\nTotal number of such non-isomorphic groups is:")
    print(number_of_groups)


solve_group_problem()