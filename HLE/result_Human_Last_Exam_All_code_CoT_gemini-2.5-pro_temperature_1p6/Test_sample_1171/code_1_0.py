def solve_rational_homotopy_vanishing():
    """
    This function determines for which k in {1, 2, ..., 9} the rational homotopy
    group pi_k(S^4 v CP^2) tensor Q vanishes.

    The logic is based on the following facts from rational homotopy theory:
    1. pi_k(A v B) tensor Q is isomorphic to (pi_k(A) tensor Q) + (pi_k(B) tensor Q)
       for simply-connected spaces A and B. This sum vanishes if and only if both
       terms vanish.
    2. The dimensions k for which pi_k(S^4) tensor Q is non-zero in the range 1-9
       are k=4 and k=7.
    3. The dimensions k for which pi_k(CP^2) tensor Q is non-zero are k=2 and k=5.
    """

    # The values of k for which pi_k(S^4) tensor Q is non-zero (for k <= 9)
    s4_nonzero_k = {4, 7}

    # The values of k for which pi_k(CP^2) tensor Q is non-zero
    cp2_nonzero_k = {2, 5}

    # We are looking for k where both rational homotopy groups are zero.
    vanishing_k = []
    for k in range(1, 10):
        # Check if the k-th rational homotopy group is non-zero for either space
        is_nonzero_for_s4 = k in s4_nonzero_k
        is_nonzero_for_cp2 = k in cp2_nonzero_k

        # The group for the wedge sum vanishes if it's zero for both components
        if not is_nonzero_for_s4 and not is_nonzero_for_cp2:
            vanishing_k.append(str(k))

    # Print the result in the specified comma-separated format
    print(",".join(vanishing_k))

solve_rational_homotopy_vanishing()