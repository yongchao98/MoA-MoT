def solve_rational_homotopy():
    """
    Finds the integers k in {1, ..., 9} for which the rational homotopy group
    pi_k(S^4 v CP^2) is trivial.
    """

    # The set of k for which the k-th rational homotopy group is non-trivial.
    # For S^4, pi_k(S^4) tensor Q is non-zero only for k=4 and k=7 in the range 1-9.
    non_zero_s4 = {4, 7}
    
    # For CP^2, pi_k(CP^2) tensor Q is non-zero only for k=2 and k=5.
    non_zero_cp2 = {2, 5}

    # For a wedge sum X = A v B, pi_k(X) tensor Q is the direct sum of the
    # rational homotopy groups of A and B. It is non-zero if either component is non-zero.
    non_zero_X = non_zero_s4.union(non_zero_cp2)

    # The set of all k we are interested in.
    k_range = set(range(1, 10))

    # The rational homotopy group vanishes for k not in the non_zero_X set.
    vanishing_k = sorted(list(k_range.difference(non_zero_X)))
    
    # Print the result in the required comma-separated format.
    print(','.join(map(str, vanishing_k)))

solve_rational_homotopy()