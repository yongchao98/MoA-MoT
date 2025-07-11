def solve_homotopy_vanishing():
    """
    Determines for which k in {1, ..., 9} the rational homotopy group
    pi_k(S^4 V CP^2) @ Q vanishes.
    """
    
    # The range of k to check
    k_range = range(1, 10)

    # For S^4, n=2. Rational homotopy is non-trivial for k = 2n and k = 4n-1.
    s4_nontrivial_k = {2 * 2, 4 * 2 - 1}
    
    # For CP^2, rational homotopy is non-trivial for k = 2 and k = 5.
    cp2_nontrivial_k = {2, 5}

    # For the wedge sum X = S^4 V CP^2, pi_k(X) @ Q is non-trivial if it is
    # non-trivial for either component.
    x_nontrivial_k = s4_nontrivial_k.union(cp2_nontrivial_k)

    # The group vanishes for all other k in the range.
    vanishing_k = [k for k in k_range if k not in x_nontrivial_k]

    # Print the result in the specified format
    print(','.join(map(str, vanishing_k)))

solve_homotopy_vanishing()