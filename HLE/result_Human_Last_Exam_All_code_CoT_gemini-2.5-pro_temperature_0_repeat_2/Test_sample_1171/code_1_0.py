def find_vanishing_rational_homotopy_groups():
    """
    Calculates for which k in {1, ..., 9} the group pi_k(S^4 v CP^2) @ Q vanishes.
    """
    # For S^4, pi_k(S^4) @ Q is non-zero for k=4 and k=7.
    s4_nonzero_k = {4, 7}

    # For CP^2, pi_k(CP^2) @ Q is non-zero for k=2 and k=5.
    cp2_nonzero_k = {2, 5}

    # For the wedge sum X = S^4 v CP^2, pi_k(X) @ Q is non-zero if k is in the union of the above sets.
    x_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)

    # We are looking for k in {1, ..., 9} for which the group vanishes.
    # These are the k's not in x_nonzero_k.
    vanishing_k = []
    for k in range(1, 10):
        if k not in x_nonzero_k:
            vanishing_k.append(str(k))

    # Print the result as a comma-separated string.
    print(",".join(vanishing_k))

find_vanishing_rational_homotopy_groups()