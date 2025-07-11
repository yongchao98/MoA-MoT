def solve_homotopy_vanishing():
    """
    This function determines for which k in {1, 2, ..., 9} the rational homotopy group
    pi_k(X) tensor Q vanishes, where X is the wedge sum of S^4 and CP^2.
    """

    # For a wedge sum of simply connected spaces like X = S^4 v CP^2,
    # the rational homotopy group pi_k(X) tensor Q vanishes if and only if
    # both pi_k(S^4) tensor Q and pi_k(CP^2) tensor Q vanish.

    # According to rational homotopy theory, the non-vanishing rational homotopy groups
    # for S^4 in the range k=1..9 are for k=4 and k=7.
    non_zero_k_for_S4 = {4, 7}

    # For CP^2, the non-vanishing rational homotopy groups are for k=2 and k=5.
    non_zero_k_for_CP2 = {2, 5}

    # The rational homotopy group of the wedge sum X is non-zero if k is in the union
    # of the above sets.
    non_zero_k_for_X = non_zero_k_for_S4.union(non_zero_k_for_CP2)

    # We are interested in k from 1 to 9.
    k_range = set(range(1, 10))

    # The group vanishes for k that are in our range but not in the non-zero set.
    vanishing_k = sorted(list(k_range.difference(non_zero_k_for_X)))

    # The final output should be a comma-separated string of these numbers.
    result_string = ",".join(map(str, vanishing_k))
    
    print(result_string)

solve_homotopy_vanishing()