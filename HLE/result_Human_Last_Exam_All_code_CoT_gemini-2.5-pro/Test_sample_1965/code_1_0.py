def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve the cube during the final 3 moves.
    This is based on the Principle of Inclusion-Exclusion.
    """
    
    # Number of possible 90-degree moves (U, U', D, D', etc.)
    N = 12

    # --- Step 1: Calculate the sizes of the individual sets ---
    
    # Let A be the set of sequences where the cube is solved after move 4.
    # To have M4*M3*M2*M1 = I, M4 must be the inverse of (M3*M2*M1).
    # We can choose M1, M2, M3 freely (N^3 ways). M4 is then fixed.
    # M5 and M6 can also be chosen freely (N^2 ways).
    # |A| = N^3 * 1 * N^2 = N^5
    count_A = N**5

    # Let B be the set of sequences where the cube is solved after move 5.
    # To have M5*M4*M3*M2*M1 = I, M5 must be the inverse of (M4*...*M1).
    # We can choose M1..M4 freely (N^4 ways). M5 is then fixed.
    # M6 can be chosen freely (N ways).
    # |B| = N^4 * 1 * N = N^5
    count_B = N**5

    # Let C be the set of sequences where the cube is solved after move 6.
    # To have M6*...*M1 = I, M6 must be the inverse of (M5*...*M1).
    # We can choose M1..M5 freely (N^5 ways). M6 is then fixed.
    # |C| = N^5 * 1 = N^5
    count_C = N**5

    # --- Step 2: Calculate the sizes of the intersections ---

    # |A intersect B|: Solved after 4 AND 5.
    # If solved after 4, the state is I. To be solved after 5, M5(I) must be I.
    # This implies M5 is the identity move, which is not in our set of 12 moves.
    # So, this is impossible.
    count_A_n_B = 0

    # |B intersect C|: Solved after 5 AND 6.
    # Similarly, if solved after 5, M6(I) must be I. This is impossible.
    count_B_n_C = 0

    # |A intersect C|: Solved after 4 AND 6.
    # Solved after 4: M4*M3*M2*M1 = I. (N^3 choices for M1-M3, M4 fixed).
    # Solved after 6: M6*M5*(State at 4) = I => M6*M5*I = I => M6*M5 = I.
    # For M6*M5=I, M6 must be the inverse of M5. We can choose M5 freely (N ways), M6 is fixed.
    # Total = (Ways for M1-M4) * (Ways for M5-M6) = (N^3 * 1) * (N * 1) = N^4
    count_A_n_C = N**4
    
    # |A intersect B intersect C|: Since A intersect B is 0, the triple intersection is also 0.
    count_A_n_B_n_C = 0

    # --- Step 3: Apply the Principle of Inclusion-Exclusion ---
    
    sum_of_singles = count_A + count_B + count_C
    sum_of_doubles = count_A_n_B + count_A_n_C + count_B_n_C
    sum_of_triples = count_A_n_B_n_C
    
    final_result = sum_of_singles - sum_of_doubles + sum_of_triples

    # --- Step 4: Print the explanation and final equation ---
    
    print("This problem can be solved using the Principle of Inclusion-Exclusion.")
    print("Let A, B, and C be the events that the cube is solved after moves 4, 5, and 6, respectively.")
    print(f"The number of possible moves at each step is {N}.\n")

    print("Calculating the size of each set:")
    print(f"|A| (solved after move 4) = {N}^5 = {count_A}")
    print(f"|B| (solved after move 5) = {N}^5 = {count_B}")
    print(f"|C| (solved after move 6) = {N}^5 = {count_C}\n")

    print("Calculating the size of the intersections:")
    print(f"|A ∩ B| (solved after 4 and 5) = {count_A_n_B} (Impossible as M5 would have to be the identity move)")
    print(f"|A ∩ C| (solved after 4 and 6) = {N}^4 = {count_A_n_C}")
    print(f"|B ∩ C| (solved after 5 and 6) = {count_B_n_C} (Impossible as M6 would have to be the identity move)")
    print(f"|A ∩ B ∩ C| (solved after 4, 5, and 6) = {count_A_n_B_n_C}\n")
    
    print("Applying the formula: |A U B U C| = |A|+|B|+|C| - (|A∩B|+|A∩C|+|B∩C|) + |A∩B∩C|")
    print("The final equation is:")
    print(f"({count_A} + {count_B} + {count_C}) - ({count_A_n_B} + {count_A_n_C} + {count_B_n_C}) + {count_A_n_B_n_C} = {final_result}")

if __name__ == "__main__":
    solve_rubiks_permutations()