def solve():
    """
    Finds the integers k in {1, ..., 9} for which the rational homotopy group
    pi_k(S^4 v CP^2) is trivial.
    """

    # The degrees k for which pi_k(S^4) tensor Q is non-trivial are 4 and 7.
    s4_nonzero_k = {4, 7}

    # The degrees k for which pi_k(CP^2) tensor Q is non-trivial are 2 and 5.
    cp2_nonzero_k = {2, 5}

    # For the wedge sum X = S^4 v CP^2, pi_k(X) tensor Q is non-trivial
    # if either of the constituent groups is non-trivial.
    # This corresponds to the union of the sets of non-trivial degrees.
    x_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)

    # We are interested in k from 1 to 9.
    k_range = set(range(1, 10))

    # The group vanishes for k where it is not non-trivial.
    # This is the set difference.
    vanishing_k = sorted(list(k_range - x_nonzero_k))

    # Print the result in the required format.
    print(','.join(map(str, vanishing_k)))

solve()