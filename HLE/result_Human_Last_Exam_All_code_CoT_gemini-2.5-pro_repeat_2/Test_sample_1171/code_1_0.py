def solve_homotopy_vanishing():
    """
    Finds the integers k in {1, ..., 9} for which the rational homotopy group
    pi_k(S^4 v CP^2) is zero.
    """
    
    # The set of k values to check.
    k_range = set(range(1, 10))
    
    # According to rational homotopy theory, for an even dimensional sphere S^n,
    # pi_k(S^n) tensor Q is non-zero only for k=n and k=2n-1.
    # For S^4, n=4, so the non-zero dimensions are k=4 and k=2*4-1=7.
    s4_nonzero_dims = {4, 7}
    
    # For the complex projective plane CP^2, the rational cohomology is Q[x]/(x^3) with |x|=2.
    # From its Sullivan minimal model, the non-zero rational homotopy groups are in degrees 2 and 5.
    cp2_nonzero_dims = {2, 5}
    
    # The rational homotopy group of the wedge sum X = S^4 v CP^2 is the direct sum
    # of the rational homotopy groups of the components.
    # pi_k(X) tensor Q = (pi_k(S^4) tensor Q) + (pi_k(CP^2) tensor Q).
    # This group is non-zero if at least one of the components is non-zero.
    x_nonzero_dims = s4_nonzero_dims.union(cp2_nonzero_dims)
    
    # The group vanishes for all other k in the range.
    vanishing_dims = sorted(list(k_range - x_nonzero_dims))
    
    # Print the result in the required comma-separated format.
    print(','.join(map(str, vanishing_dims)))

solve_homotopy_vanishing()