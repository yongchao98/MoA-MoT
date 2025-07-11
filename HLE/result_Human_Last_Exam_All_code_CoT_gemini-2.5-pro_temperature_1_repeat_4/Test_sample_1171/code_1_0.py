def solve_homotopy_vanishing():
    """
    Calculates for which k in {1, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) tensor Q vanishes.
    """

    # The full range of k values
    k_range = set(range(1, 10))

    # For S^n (n even, n=2m), pi_k(S^n) tensor Q is non-zero for k=n and k=2n-1.
    # For S^4, n=4, m=2. Non-zero for k=4 and k=2*4-1=7.
    s4_nonzero_k = {4, 7}

    # For CP^n, pi_k(CP^n) tensor Q is non-zero for k=2 and k=2n+1.
    # For CP^2, n=2. Non-zero for k=2 and k=2*2+1=5.
    cp2_nonzero_k = {2, 5}

    # For the wedge sum X = S^4 v CP^2, pi_k(X) tensor Q is non-zero
    # if it's non-zero for either component.
    x_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)

    # The group vanishes for all other k in the range.
    vanishing_k = sorted(list(k_range.difference(x_nonzero_k)))

    # Print the result in the specified format.
    print(",".join(map(str, vanishing_k)))

solve_homotopy_vanishing()