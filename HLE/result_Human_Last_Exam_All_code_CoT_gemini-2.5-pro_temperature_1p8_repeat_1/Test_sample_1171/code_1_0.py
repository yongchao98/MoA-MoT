def solve_rational_homotopy_vanishing():
    """
    Calculates for which k in {1, ..., 9} the group pi_k(S^4 v CP^2) @ Q vanishes.
    """
    
    # The group pi_k(S^4 v CP^2) @ Q vanishes if and only if
    # both pi_k(S^4) @ Q and pi_k(CP^2) @ Q vanish.
    # We first find the k for which these groups are NON-ZERO.

    # For an even sphere S^{2m}, non-zero rational homotopy is at k = 2m and k = 4m-1.
    # For S^4, m=2.
    s4_m = 2
    s4_nonzero_k = {2 * s4_m, 4 * s4_m - 1}
    print(f"k for which pi_k(S^4) @ Q is non-zero: {sorted(list(s4_nonzero_k))}")

    # For CP^n, non-zero rational homotopy is at k = 2 and k = 2n+1.
    # For CP^2, n=2.
    cp2_n = 2
    cp2_nonzero_k = {2, 2 * cp2_n + 1}
    print(f"k for which pi_k(CP^2) @ Q is non-zero: {sorted(list(cp2_nonzero_k))}")

    # The group for the wedge sum X is non-zero if k is in the union of the above sets.
    total_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)
    print(f"k for which pi_k(X) @ Q is non-zero: {sorted(list(total_nonzero_k))}")

    # We want k where the group vanishes. This is the complement in {1, ..., 9}.
    k_range = set(range(1, 10))
    vanishing_k = k_range.difference(total_nonzero_k)
    
    # Sort and format the final answer.
    result_string = ",".join(map(str, sorted(list(vanishing_k))))
    
    print("\nFinal Answer:")
    print(result_string)

solve_rational_homotopy_vanishing()