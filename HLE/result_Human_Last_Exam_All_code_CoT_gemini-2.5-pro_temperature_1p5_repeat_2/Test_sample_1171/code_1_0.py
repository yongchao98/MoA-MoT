def solve_homotopy():
    """
    Finds the values of k in {1, 2, ..., 9} for which the rational
    homotopy group pi_k(S^4 v CP^2) vanishes.
    """
    
    # The set of k values to check.
    k_range = range(1, 10)
    
    # Determine the values of k for which pi_k(S^4) @ Q is non-zero.
    # For S^{2m}, non-zero rational homotopy groups are in degrees 2m and 4m-1.
    # For S^4, m=2. So non-zero degrees are k=4 and k=7.
    s4_m = 2
    s4_nonzero_k = {2 * s4_m, 4 * s4_m - 1}
    
    # Determine the values of k for which pi_k(CP^2) @ Q is non-zero.
    # For CP^n, non-zero rational homotopy groups are in degrees 2 and 2n+1.
    # For CP^2, n=2. So non-zero degrees are k=2 and k=5.
    cp2_n = 2
    cp2_nonzero_k = {2, 2 * cp2_n + 1}
    
    # pi_k(S^4 v CP^2) @ Q is non-zero if either of the component groups is non-zero.
    total_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)
    
    # We want the values of k for which the group vanishes.
    vanishing_k = []
    for k in k_range:
        if k not in total_nonzero_k:
            vanishing_k.append(k)
            
    # Print the final list of k values, comma-separated.
    print(",".join(map(str, vanishing_k)))

solve_homotopy()