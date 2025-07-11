def solve_homotopy_vanishing():
    """
    Determines for which k in {1, 2, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) vanishes.

    The rational homotopy group of a wedge sum of simply connected spaces is the
    direct sum of the rational homotopy groups of the components.
    pi_k(X) @ Q = (pi_k(S^4) @ Q) + (pi_k(CP^2) @ Q).
    This group vanishes if and only if both summands vanish.

    1. Rational homotopy groups of S^{2n}:
       pi_k(S^{2n}) @ Q is non-zero only for k = 2n and k = 4n - 1.
       For S^4, n=2. Non-zero for k=4 and k=7.

    2. Rational homotopy groups of CP^n:
       pi_k(CP^n) @ Q is non-zero only for k=2 and k=2n+1.
       For CP^2, n=2. Non-zero for k=2 and k=5.

    We need to find k in {1, ..., 9} such that k is not in {4, 7} and k is not in {2, 5}.
    """

    # k values for which pi_k(S^4) @ Q is non-zero
    s4_nonzero_k = {4, 7}

    # k values for which pi_k(CP^2) @ Q is non-zero
    cp2_nonzero_k = {2, 5}

    # The union of these sets gives all k for which pi_k(X) @ Q is non-zero
    total_nonzero_k = s4_nonzero_k.union(cp2_nonzero_k)

    vanishing_k = []
    for k in range(1, 10):
        if k not in total_nonzero_k:
            vanishing_k.append(str(k))

    print(",".join(vanishing_k))

solve_homotopy_vanishing()