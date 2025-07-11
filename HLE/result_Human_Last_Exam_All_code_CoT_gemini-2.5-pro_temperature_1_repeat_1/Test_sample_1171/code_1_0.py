def solve_homotopy_vanishing():
    """
    Finds for which k in {1, ..., 9} the group pi_k(S^4 v CP^2) @ Q vanishes.
    """
    # The set of k values we are interested in.
    k_values = set(range(1, 10))

    # According to rational homotopy theory, the non-trivial rational homotopy groups of S^n (for n even)
    # occur at dimensions n and 2n-1. For S^4, n=4.
    n_s4 = 4
    non_trivial_k_s4 = {n_s4, 2 * n_s4 - 1}

    # The non-trivial rational homotopy groups of CP^n occur at dimensions 2 and 2n+1.
    # For CP^2, n=2.
    n_cp2 = 2
    non_trivial_k_cp2 = {2, 2 * n_cp2 + 1}

    # For the wedge sum X = S^4 v CP^2, pi_k(X) @ Q is non-trivial if pi_k(S^4) @ Q or pi_k(CP^2) @ Q is non-trivial.
    # So, the set of non-trivial k is the union of the two sets above.
    non_trivial_k_X = non_trivial_k_s4.union(non_trivial_k_cp2)

    # The group vanishes for all other k in the specified range.
    # This is the set difference between the total set of k and the non-trivial set.
    vanishing_k = sorted(list(k_values.difference(non_trivial_k_X)))

    # Print the result in the specified format.
    print(','.join(map(str, vanishing_k)))

solve_homotopy_vanishing()