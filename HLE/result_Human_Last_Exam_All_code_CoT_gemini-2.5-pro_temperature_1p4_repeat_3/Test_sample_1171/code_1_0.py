def solve_homotopy_vanishing():
    """
    Finds for which k in {1,...,9} the rational homotopy group
    pi_k(S^4 v CP^2) tensor Q vanishes.
    """

    # The range of k we are interested in.
    k_range = set(range(1, 10))

    # According to rational homotopy theory, for a wedge sum of simply-connected
    # spaces X = A v B, we have pi_k(X) @ Q = (pi_k(A) @ Q) + (pi_k(B) @ Q).
    # This vanishes if and only if both summands vanish.
    # We first identify the dimensions 'k' where the groups are NON-trivial.

    # For the 4-sphere S^4, the non-trivial rational homotopy groups are
    # for k=4 and k=4*2-1=7.
    nontrivial_k_s4 = {4, 7}

    # For the complex projective plane CP^2, the non-trivial rational homotopy
    # groups are for k=2 and k=5.
    nontrivial_k_cp2 = {2, 5}

    # The rational homotopy group of the wedge sum X is non-trivial if either
    # component's group is non-trivial.
    nontrivial_k_X = nontrivial_k_s4.union(nontrivial_k_cp2)

    # The group pi_k(X) @ Q vanishes for all other k in the specified range.
    vanishing_k = sorted(list(k_range.difference(nontrivial_k_X)))

    # Print the result in the required format.
    print(','.join(map(str, vanishing_k)))

solve_homotopy_vanishing()