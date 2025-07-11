def find_vanishing_rational_homotopy_groups():
    """
    This function calculates for which k in {1, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) @ Q vanishes.

    The rational homotopy group of a wedge sum is the direct sum of the
    rational homotopy groups of its components.
    pi_k(X) @ Q = (pi_k(S^4) @ Q) + (pi_k(CP^2) @ Q)
    This group vanishes if and only if both summands vanish.

    We first identify the k for which the groups are non-zero.
    """

    # For S^4, pi_k(S^4) @ Q is non-zero for k=4 and k=2*4-1=7.
    s4_non_vanishing_k = {4, 7}
    print(f"For S^4, the rational homotopy groups pi_k(S^4) @ Q are non-zero for k = {sorted(list(s4_non_vanishing_k))}")

    # For CP^2, pi_k(CP^2) @ Q is non-zero for k=2 and k=2*2+1=5.
    cp2_non_vanishing_k = {2, 5}
    print(f"For CP^2, the rational homotopy groups pi_k(CP^2) @ Q are non-zero for k = {sorted(list(cp2_non_vanishing_k))}")

    # For the wedge sum X = S^4 v CP^2, pi_k(X) @ Q is non-zero if k is in the union of the above sets.
    total_non_vanishing_k = s4_non_vanishing_k.union(cp2_non_vanishing_k)
    print(f"Therefore, for X = S^4 v CP^2, pi_k(X) @ Q is non-zero for k in {sorted(list(total_non_vanishing_k))}")

    # We are looking for k in {1, ..., 9} for which the group vanishes.
    # These are the k's not in the non-vanishing set.
    vanishing_k = []
    for k in range(1, 10):
        if k not in total_non_vanishing_k:
            vanishing_k.append(k)

    print("\nThus, pi_k(X) @ Q vanishes for the following values of k in {1, 2, ..., 9}:")
    # The final result is printed in the required comma-separated format.
    print(','.join(map(str, vanishing_k)))

find_vanishing_rational_homotopy_groups()