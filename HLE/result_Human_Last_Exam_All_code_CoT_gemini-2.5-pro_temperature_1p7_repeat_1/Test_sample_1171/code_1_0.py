def solve_homotopy_problem():
    """
    This function determines for which k in {1, ..., 9} the rational homotopy
    group pi_k(S^4 v CP^2) vanishes.
    """

    # The rational homotopy groups of the wedge sum are the direct sum of the
    # rational homotopy groups of the components, since S^4 and CP^2 are
    # simply connected.
    # pi_k(X) (x) Q = (pi_k(S^4) (x) Q) + (pi_k(CP^2) (x) Q)
    # This group vanishes if and only if both summands vanish.

    # According to rational homotopy theory, pi_k(S^n) (x) Q is non-zero only for
    # k=n and k=2n-1 for n even. For S^4 (n=4), k = 4 and k = 2*4-1 = 7.
    s4_nontrivial_k = {4, 7}

    # pi_k(CP^n) (x) Q is non-zero only for k=2 and k=2n+1.
    # For CP^2 (n=2), k = 2 and k = 2*2+1 = 5.
    cp2_nontrivial_k = {2, 5}

    # The rational homotopy group of the wedge sum X is non-trivial if either
    # component's group is non-trivial.
    x_nontrivial_k = s4_nontrivial_k.union(cp2_nontrivial_k)

    # We are interested in k from 1 to 9.
    k_range = set(range(1, 10))

    # The k for which the group vanishes are those not in the non-trivial set.
    vanishing_k = sorted(list(k_range - x_nontrivial_k))

    # Print the result in the specified format
    print(",".join(map(str, vanishing_k)))

solve_homotopy_problem()