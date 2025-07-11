def solve_homotopy_problem():
    """
    Determines for which k in {1, 2, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) tensor Q vanishes.

    According to the rational Hilton-Milnor theorem, for simply connected spaces A and B,
    pi_k(A v B) tensor Q is isomorphic to (pi_k(A) tensor Q) + (pi_k(B) tensor Q).
    This group vanishes if and only if both component groups vanish.
    """

    # Non-zero rational homotopy groups for S^4 occur at k = 4, 7.
    # For a sphere S^(2m), they are at k=2m and k=4m-1. Here m=2.
    s4_nonzero_k = {4, 7}

    # Non-zero rational homotopy groups for CP^2 occur at k = 2, 5.
    # For CP^n, they are at k=2 and k=2n+1. Here n=2.
    cp2_nonzero_k = {2, 5}

    # The group for the wedge sum is non-zero if k is in the union of the above sets.
    total_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)

    # We want the values of k for which the group vanishes.
    vanishing_k = []
    for k in range(1, 10):
        if k not in total_nonzero_k:
            vanishing_k.append(k)

    # Print the final result as a comma-separated string.
    print(",".join(map(str, vanishing_k)))

solve_homotopy_problem()