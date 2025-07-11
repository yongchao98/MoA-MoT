def solve_rational_homotopy():
    """
    Calculates for which k in {1, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) @ Q vanishes.
    """

    # Degrees k for which pi_k(S^4) @ Q is non-trivial (for k <= 9)
    # For an even sphere S^{2n}, non-trivial rational homotopy is in degrees 2n and 4n-1.
    # For S^4, n=2, so degrees are 4 and 7.
    s4_nonzero_k = {4, 7}

    # Degrees k for which pi_k(CP^2) @ Q is non-trivial (for k <= 9)
    # For CP^n, non-trivial rational homotopy is in degrees 2 and 2n+1.
    # For CP^2, n=2, so degrees are 2 and 5.
    cp2_nonzero_k = {2, 5}

    # For the wedge sum X = S^4 v CP^2, pi_k(X) @ Q is the direct sum.
    # It is non-trivial if either summand is non-trivial.
    # So, the set of non-trivial degrees for X is the union of the sets above.
    x_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)

    # We are looking for k in {1, ..., 9} where the group *vanishes*.
    # This is the complement of x_nonzero_k in the set {1, ..., 9}.
    all_k = set(range(1, 10))
    vanishing_k = sorted(list(all_k - x_nonzero_k))

    # Print the result in the required format
    print(','.join(map(str, vanishing_k)))

solve_rational_homotopy()