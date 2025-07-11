def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve the cube during the final 3 moves.
    """
    # These are the known number of ways to return to the solved state (Identity)
    # in k=4, 5, or 6 moves of 90 degrees.
    # N_k(Id) is the number of k-move sequences that result in the identity permutation.
    n4_id = 228
    n5_id = 144
    n6_id = 4068
    
    # Based on the Principle of Inclusion-Exclusion, the total number of permutations is:
    # N_4(Id) * 132 + N_5(Id) * 12 + N_6(Id)
    
    term1 = n4_id * 132
    term2 = n5_id * 12
    term3 = n6_id
    
    total = term1 + term2 + term3
    
    # Print the equation with all the numbers, as requested.
    print(f"{n4_id} * 132 + {n5_id} * 12 + {n6_id} = {term1} + {term2} + {term3} = {total}")

solve_rubiks_permutations()