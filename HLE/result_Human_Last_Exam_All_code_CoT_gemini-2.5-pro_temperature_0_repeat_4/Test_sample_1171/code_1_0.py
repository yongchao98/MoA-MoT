def solve():
    """
    Determines for which k in {1, ..., 9} the group pi_k(S^4 v CP^2) tensor Q vanishes.
    """
    
    # Set of k where the rational homotopy groups are non-zero for S^4
    # For S^n with n even, non-zero for k=n and k=2n-1. Here n=4.
    non_vanishing_s4 = {4, 2*4 - 1} # {4, 7}

    # Set of k where the rational homotopy groups are non-zero for CP^2
    # For CP^n, non-zero for k=2 and k=2n+1. Here n=2.
    non_vanishing_cp2 = {2, 2*2 + 1} # {2, 5}

    # The rational homotopy group of the wedge sum is the direct sum of the components.
    # It vanishes if and only if both components' rational homotopy groups vanish.
    # So, the non-vanishing k for the wedge sum is the union of the individual sets.
    non_vanishing_x = non_vanishing_s4.union(non_vanishing_cp2)

    vanishing_k = []
    for k in range(1, 10):
        if k not in non_vanishing_x:
            vanishing_k.append(k)
            
    # Print the result in the required format
    print(",".join(map(str, vanishing_k)))

solve()