def solve_homotopy_vanishing():
    """
    This function calculates the integers k in {1, ..., 9} for which
    the rational homotopy group pi_k(S^4 v CP^2) @ Q vanishes.
    """

    # The set of k values to check.
    k_range = set(range(1, 10))

    # According to rational homotopy theory, the degrees of the generators
    # of the minimal Sullivan model for X = A v B are the union of the
    # degrees of the generators for A and B.

    # For S^4, the generators are in degrees 4 and 7.
    non_vanishing_s4 = {4, 7}

    # For CP^2, the generators are in degrees 2 and 5.
    non_vanishing_cp2 = {2, 5}

    # For the wedge sum X = S^4 v CP^2, the non-vanishing rational homotopy groups
    # correspond to the union of these sets of degrees.
    non_vanishing_x = non_vanishing_s4.union(non_vanishing_cp2)

    # The set of k for which the group vanishes is the complement
    # of the non-vanishing set in the given range.
    vanishing_k = sorted(list(k_range - non_vanishing_x))

    # Print the result in the specified comma-separated format.
    print(",".join(map(str, vanishing_k)))

solve_homotopy_vanishing()