def solve_rational_homotopy_vanishing():
    """
    Calculates for which k in {1, ..., 9} the group pi_k(S^4 v CP^2) @ Q vanishes.
    
    The rational homotopy group of a wedge sum of simply-connected spaces is the
    direct sum of the individual rational homotopy groups. Therefore, the group
    vanishes if and only if it vanishes for both constituent spaces.

    1. For S^4, the non-zero rational homotopy groups are at k=4 and k=7.
    2. For CP^2, the non-zero rational homotopy groups are at k=2 and k=5.
    
    We seek k in {1,...,9} which are in neither of these sets.
    """
    k_range = set(range(1, 10))

    # k for which pi_k(S^4) @ Q is non-zero
    s4_nonzero = {4, 7}
    
    # k for which pi_k(CP^2) @ Q is non-zero
    cp2_nonzero = {2, 5}
    
    # k for which the group for S^4 vanishes
    s4_zero = k_range.difference(s4_nonzero)
    
    # k for which the group for CP^2 vanishes
    cp2_zero = k_range.difference(cp2_nonzero)
    
    # The result is the intersection of these two sets
    result_k = sorted(list(s4_zero.intersection(cp2_zero)))
    
    # Print the result in the specified format
    print(','.join(map(str, result_k)))

solve_rational_homotopy_vanishing()