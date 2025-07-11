def find_vanishing_homotopy_groups():
    """
    Finds the values of k in {1, 2, ..., 9} for which the rational homotopy group
    pi_k(S^4 v CP^2) tensor Q vanishes.

    The rational homotopy group of the wedge sum vanishes if and only if the
    rational homotopy groups of both constituent spaces vanish.

    We first identify the k for which the groups are non-zero.
    """

    # For S^4, pi_k(S^4) tensor Q is non-zero for k = 4 and k = 2*4 - 1 = 7.
    nonzero_k_for_s4 = {4, 7}

    # For CP^2, pi_k(CP^2) tensor Q is non-zero for k = 2 and k = 5.
    nonzero_k_for_cp2 = {2, 5}

    # The group for the wedge sum X is non-zero if k is in the union of the above sets.
    nonzero_k_for_X = nonzero_k_for_s4.union(nonzero_k_for_cp2)

    # The set of k we are interested in is {1, 2, ..., 9}.
    k_range = set(range(1, 10))

    # The group vanishes for k in the complement of the non-zero set.
    vanishing_k = sorted(list(k_range - nonzero_k_for_X))

    # Print the result as a comma-separated string.
    print(','.join(map(str, vanishing_k)))

find_vanishing_homotopy_groups()