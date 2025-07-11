def find_vanishing_rational_homotopy_groups():
    """
    Determines for which k in {1, 2, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) tensor Q vanishes.

    The rational homotopy group pi_k(S^4 v CP^2) tensor Q is the direct sum
    of (pi_k(S^4) tensor Q) and (pi_k(CP^2) tensor Q). It vanishes if and
    only if both summands vanish.
    """

    # For S^4, the non-vanishing rational homotopy groups are for k=4 and k=7.
    s4_nonzero_k = {4, 7}

    # For CP^2, the non-vanishing rational homotopy groups are for k=2 and k=5.
    cp2_nonzero_k = {2, 5}

    vanishing_k_values = []

    # We check for k from 1 to 9.
    for k in range(1, 10):
        # Check if the rational homotopy group for S^4 vanishes at k.
        pi_k_s4_is_zero = k not in s4_nonzero_k
        
        # Check if the rational homotopy group for CP^2 vanishes at k.
        pi_k_cp2_is_zero = k not in cp2_nonzero_k
        
        # The group for the wedge sum vanishes if and only if both components vanish.
        if pi_k_s4_is_zero and pi_k_cp2_is_zero:
            vanishing_k_values.append(str(k))

    # Print the result in the required format.
    print(",".join(vanishing_k_values))

find_vanishing_rational_homotopy_groups()