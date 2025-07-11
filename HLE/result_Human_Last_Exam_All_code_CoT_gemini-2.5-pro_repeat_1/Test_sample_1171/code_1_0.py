def solve_homotopy_vanishing():
    """
    Calculates for which k in {1, ..., 9} the group pi_k(S^4 v CP^2) tensor Q vanishes.
    """
    # The set of k values to consider
    k_values = set(range(1, 10))

    # The dimensions 'k' for which the rational homotopy groups of S^4 are non-zero.
    # For S^4, n=2, so non-vanishing at k=2n=4 and k=4n-1=7.
    s4_non_vanishing = {4, 7}

    # The dimensions 'k' for which the rational homotopy groups of CP^2 are non-zero.
    # For CP^2, n=2, so non-vanishing at k=2 and k=2n+1=5.
    cp2_non_vanishing = {2, 5}

    # For the wedge sum X = S^4 v CP^2, pi_k(X) tensor Q is non-zero if either
    # of the component groups is non-zero. The set of non-vanishing dimensions for X is the union.
    x_non_vanishing = s4_non_vanishing.union(cp2_non_vanishing)

    # We want the k values for which the group vanishes.
    # This is the set difference between all k_values and the non-vanishing ones.
    vanishing_k = sorted(list(k_values.difference(x_non_vanishing)))

    # Format the result as a comma-separated string.
    result_string = ",".join(map(str, vanishing_k))
    print(result_string)

solve_homotopy_vanishing()