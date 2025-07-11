def solve():
    """
    Finds the integers k in {1, 2, ..., 9} for which the rational homotopy group
    pi_k(S^4 v CP^2) tensor Q vanishes.
    """
    
    # According to rational homotopy theory, the non-vanishing rational homotopy groups
    # for S^4 and CP^2 are at specific dimensions.
    # For S^n (n even), pi_k(S^n) tensor Q is non-zero for k=n and k=2n-1.
    # For S^4, this is k=4 and k=7.
    s4_nonzero_k = {4, 7}
    
    # For CP^n, pi_k(CP^n) tensor Q is non-zero for k=2 and k=2n+1.
    # For CP^2, this is k=2 and k=5.
    cp2_nonzero_k = {2, 5}
    
    # For k>=2, pi_k(S^4 v CP^2) tensor Q is the direct sum of the individual groups.
    # It is non-zero if either of the individual groups is non-zero.
    # The set of k where the group is non-zero is the union of the two sets above.
    total_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)
    
    vanishing_k = []
    # We check for k from 1 to 9.
    for k in range(1, 10):
        # For k=1, pi_1(S^4 v CP^2) is trivial because both spaces are simply connected.
        # For k>=2, the group vanishes if k is not in the set of non-zero dimensions.
        if k == 1 or k not in total_nonzero_k:
            vanishing_k.append(k)
            
    # Print the result in the specified format
    print(','.join(map(str, vanishing_k)))

solve()