def solve():
    """
    Finds for which k in {1, ..., 9} the group pi_k(S^4 v CP^2) tensor Q vanishes.
    """
    
    # For S^4 (an even sphere S^{2n} with n=2), the rational homotopy groups
    # are non-zero only for k=2n and k=4n-1.
    non_vanishing_s4 = {2*2, 4*2 - 1}  # {4, 7}
    
    # For CP^2, the rational homotopy groups are non-zero for k=2 (by Hurewicz)
    # and k=5 (since pi_k(CP^2) is isomorphic to pi_k(S^5) for k>2).
    non_vanishing_cp2 = {2, 5}
    
    # The rational homotopy group of the wedge sum is the direct sum of the
    # individual rational homotopy groups. It vanishes if and only if both parts vanish.
    # So, the group is non-zero if k is in either of the non-vanishing sets.
    non_vanishing_x = non_vanishing_s4.union(non_vanishing_cp2)
    
    vanishing_k = []
    for k in range(1, 10):
        if k not in non_vanishing_x:
            vanishing_k.append(k)
            
    # Print the result in the required comma-separated format.
    print(','.join(map(str, vanishing_k)))

solve()