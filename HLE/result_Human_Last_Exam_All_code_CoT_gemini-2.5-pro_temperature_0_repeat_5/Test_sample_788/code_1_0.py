def solve():
    """
    Calculates the number of equivalence classes for the peg game.

    The total number of classes is the product of two invariants:
    1. Parity Invariant: Splits configurations into "even" and "odd". (2 classes)
    2. Coordinate Vector Invariant: The number of cosets of a subspace W in Z_2^4.
    """

    # The number of classes from the parity invariant is 2.
    num_parity_classes = 2

    # --- Calculate the number of classes from the second invariant ---

    # The vector space is Z_2^4. Addition is XOR.
    def add_vectors(v1, v2):
        return tuple((v1[i] ^ v2[i]) for i in range(4))

    # The vectors in Z_2^4 corresponding to the basic horizontal and vertical moves.
    # These are derived from the definition of the coordinate parity vector.
    phi_mh = (1, 0, 1, 1)  # Vector for a horizontal move at the origin
    phi_mv = (0, 1, 1, 1)  # Vector for a vertical move at the origin
    zero_vector = (0, 0, 0, 0)

    # Generate the subspace W spanned by the move vectors.
    # The elements are 0, phi_mh, phi_mv, and phi_mh + phi_mv.
    W = {
        zero_vector,
        phi_mh,
        phi_mv,
        add_vectors(phi_mh, phi_mv)
    }
    size_W = len(W)

    # The total space is Z_2^4, which has 2^4 elements.
    size_Z2_4 = 2**4

    # The number of cosets is the size of the total space divided by the size of the subspace.
    num_coset_classes = size_Z2_4 // size_W

    # The total number of equivalence classes is the product of the two.
    total_classes = num_parity_classes * num_coset_classes

    # Print the final equation as requested.
    print(f"The total number of equivalence classes is the product of the number of parity classes ({num_parity_classes}) and the number of coordinate vector classes ({num_coset_classes}).")
    print(f"{num_parity_classes} * {num_coset_classes} = {total_classes}")

solve()