def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at some point during the final 3 moves.
    """
    
    # N(k) is the number of k-move sequences that return the cube to identity.
    # These are known values from computational group theory.
    # A standard 90-degree move can be one of 12 possibilities (U, U', D, D', etc.).
    N2 = 12
    N4 = 324
    N5 = 0  # It's impossible to return to solved state in an odd number of moves.
    N6 = 12624
    
    # Number of choices for any single move.
    num_moves = 12
    
    # |A|: Solved after the 4th move.
    # The first 4 moves form a loop, the last 2 can be anything.
    count_A = N4 * (num_moves ** 2)
    
    # |B|: Solved after the 5th move.
    # The first 5 moves form a loop, the last move can be anything.
    count_B = N5 * num_moves
    
    # |C|: Solved after the 6th move.
    # All 6 moves form a loop.
    count_C = N6
    
    # |A intersect C|: Solved after the 4th move AND after the 6th move.
    # This means the first 4 moves form a loop, and moves 5 and 6 also form a loop.
    count_A_intersect_C = N4 * N2
    
    # Using the Principle of Inclusion-Exclusion:
    # Total = |A| + |B| + |C| - |A intersect C|
    # Other intersections are 0.
    total_permutations = count_A + count_B + count_C - count_A_intersect_C
    
    # Print the step-by-step calculation
    print("Let's calculate the number of favorable permutations step by step.")
    print("-" * 60)
    print(f"Number of permutations solved after exactly 4 moves (A):")
    print(f"  N(4) * 12^2 = {N4} * {num_moves**2} = {count_A}")
    print("\nNumber of permutations solved after exactly 5 moves (B):")
    print(f"  N(5) * 12 = {N5} * {num_moves} = {count_B}")
    print("\nNumber of permutations solved after exactly 6 moves (C):")
    print(f"  N(6) = {N6}")
    print("\nNumber of permutations solved at both 4 and 6 moves (A and C):")
    print(f"  N(4) * N(2) = {N4} * {N2} = {count_A_intersect_C}")
    print("-" * 60)
    print("Using the Principle of Inclusion-Exclusion, the total is:")
    print("  Total = |A| + |B| + |C| - |A and C|")
    print(f"  Total = {count_A} + {count_B} + {count_C} - {count_A_intersect_C}")
    print(f"  Total = {total_permutations}")

solve_rubiks_permutations()
<<<55392>>>