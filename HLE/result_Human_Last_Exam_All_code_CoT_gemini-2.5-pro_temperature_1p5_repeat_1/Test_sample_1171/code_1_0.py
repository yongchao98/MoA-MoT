def solve_homotopy_vanishing():
    """
    Finds for which k in {1, ..., 9} the group pi_k(S^4 v CP^2) @ Q vanishes.
    
    The rational homotopy group pi_k(A v B) @ Q is isomorphic to
    (pi_k(A) @ Q) + (pi_k(B) @ Q) for simply connected A, B.
    This vanishes if and only if both components vanish.
    
    1. For S^n, pi_k(S^n) @ Q is non-zero only for:
       - k = n
       - k = 2n-1, if n is even.
       For S^4 (n=4, even), non-zero for k=4 and k=2*4-1=7.

    2. For CP^n, pi_k(CP^n) @ Q is non-zero only for:
       - k = 2
       - k = 2n+1
       For CP^2 (n=2), non-zero for k=2 and k=2*2+1=5.
    
    So, pi_k(S^4 v CP^2) @ Q is non-zero if k is in {4, 7} U {2, 5} = {2, 4, 5, 7}.
    We want the k values for which the group *vanishes*.
    """
    
    k_range = range(1, 10)
    
    # k values for which pi_k(S^4) @ Q is non-zero
    n_sphere = 4
    s4_nonzero_k = {n_sphere}
    if n_sphere % 2 == 0:
        s4_nonzero_k.add(2 * n_sphere - 1)
        
    # k values for which pi_k(CP^2) @ Q is non-zero
    n_cp = 2
    cp2_nonzero_k = {2, 2 * n_cp + 1}
    
    # The total set of k for which pi_k(X) @ Q is non-zero
    x_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)
    
    vanishing_k = []
    for k in k_range:
        if k not in x_nonzero_k:
            vanishing_k.append(k)
            
    # Print the result in the specified comma-separated format
    print(','.join(map(str, vanishing_k)))

solve_homotopy_vanishing()