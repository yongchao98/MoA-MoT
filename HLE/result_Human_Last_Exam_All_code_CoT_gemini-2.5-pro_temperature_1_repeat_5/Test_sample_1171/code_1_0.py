def solve_homotopy_vanishing():
    """
    Finds the integers k in {1, ..., 9} for which pi_k(S^4 v CP^2) tensor Q vanishes.
    """
    # According to rational homotopy theory, pi_k(S^4 v CP^2) tensor Q is non-zero
    # if and only if either pi_k(S^4) tensor Q or pi_k(CP^2) tensor Q is non-zero.

    # Non-vanishing rational homotopy groups for S^4: k = 4, 2*4-1 = 7.
    s4_nonzero_k = {4, 7}

    # Non-vanishing rational homotopy groups for CP^2: k = 2, 2*2+1 = 5.
    cp2_nonzero_k = {2, 5}

    # The union of these sets gives the k for which pi_k(X) tensor Q is non-zero.
    x_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)

    # The range of k we are interested in.
    k_range = set(range(1, 10))

    # The set of k for which the group vanishes is the difference.
    vanishing_k = k_range.difference(x_nonzero_k)

    # Sort the results and format for printing.
    result = sorted(list(vanishing_k))
    print(",".join(map(str, result)))

solve_homotopy_vanishing()