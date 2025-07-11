def solve_rational_homotopy_vanishing():
    """
    Calculates for which k in {1, ..., 9} the group pi_k(S^4 v CP^2) tensor Q vanishes.
    """
    k_range = set(range(1, 10))

    # For S^4, pi_k(S^4) tensor Q is non-zero for k=4 and k=2*4-1=7.
    s4_nonzero_k = {4, 7}

    # For CP^2, pi_k(CP^2) tensor Q is non-zero for k=2 and k=2*2+1=5.
    cp2_nonzero_k = {2, 5}

    # For the wedge sum X = S^4 v CP^2, pi_k(X) tensor Q is non-zero if it is
    # non-zero for either component.
    x_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)

    # The group vanishes for k values in our range that are not in the non-zero set.
    vanishing_k = sorted(list(k_range - x_nonzero_k))

    # Print the final result as a comma-separated string.
    print(",".join(map(str, vanishing_k)))

solve_rational_homotopy_vanishing()