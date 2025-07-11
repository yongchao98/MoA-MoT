import math

def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that result in a solved cube
    at some point during the final 3 moves, using the Principle of Inclusion-Exclusion.
    """
    
    # Number of possible 90-degree moves (6 faces * 2 directions)
    num_moves = 12
    
    # |A|: Number of sequences solved after 4 moves.
    # The first 3 moves are free, the 4th is fixed. The last 2 are free.
    # |A| = 12^3 * 1 * 12^2 = 12^5
    A = num_moves**5
    
    # |B|: Number of sequences solved after 5 moves.
    # The first 4 moves are free, the 5th is fixed. The 6th is free.
    # |B| = 12^4 * 1 * 12^1 = 12^5
    B = num_moves**5
    
    # |C|: Number of sequences solved after 6 moves.
    # The first 5 moves are free, the 6th is fixed.
    # |C| = 12^5 * 1 = 12^5
    C = num_moves**5
    
    # |A intersect B|: Solved at move 4 and 5. Impossible, as a single 90-degree turn
    # cannot result in a solved cube from a solved state.
    A_int_B = 0
    
    # |B intersect C|: Solved at move 5 and 6. Impossible for the same reason.
    B_int_C = 0
    
    # |A intersect C|: Solved at move 4 and 6. This requires the first 4 moves to be an identity sequence,
    # and moves 5 and 6 to be inverses of each other.
    # Count = 12^3 (M1-M3) * 1 (M4) * 12 (M5) * 1 (M6) = 12^4
    A_int_C = num_moves**4
    
    # |A intersect B intersect C|: Must be 0 as it's a subset of |A intersect B|.
    A_int_B_int_C = 0
    
    # Apply the Principle of Inclusion-Exclusion
    total = A + B + C - (A_int_B + A_int_C + B_int_C) + A_int_B_int_C
    
    # Print the final equation with each component calculated.
    print("The total number of favorable permutations is calculated as:")
    print(f"|A| + |B| + |C| - (|A ∩ B| + |A ∩ C| + |B ∩ C|) + |A ∩ B ∩ C|")
    print(f"= {A} + {B} + {C} - ({A_int_B} + {A_int_C} + {B_int_C}) + {A_int_B_int_C}")
    print(f"= {total}")

solve_rubiks_permutations()
<<<725760>>>