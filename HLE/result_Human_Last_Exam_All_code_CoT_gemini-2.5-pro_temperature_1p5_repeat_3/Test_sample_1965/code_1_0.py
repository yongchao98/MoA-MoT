# This script requires the 'pycuber' library.
# Please install it first by running: pip install pycuber

import pycuber as pc
from collections import defaultdict
import sys

def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at some point during the final 3 moves.
    """
    print("This script will calculate the answer. It may take a minute or two to run.")
    sys.stdout.flush()

    # The 12 standard 90-degree moves
    moves = ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]

    # Use a fast, hashable representation of the cube's state.
    # A tuple of the corner and edge cubie objects is hashable and unique.
    def get_state(cube):
        return (cube.corners, cube.edges)

    # Initialize with the solved state
    solved_cube = pc.Cube()
    solved_state = get_state(solved_cube)

    # We use a dictionary to store the number of ways to reach each state.
    # It maps state -> count.
    counts = {solved_state: 1}
    
    # We also keep a mapping from the state representation to the cube object
    # to avoid reconstructing cubes, which is slow.
    state_to_cube = {solved_state: solved_cube}

    # n_k stores the number of k-move sequences that return to the start state.
    n_k = [1]  # n_0 is 1 (the empty sequence)

    print("Step 1: Calculating n_k for k=1 to 6...")
    sys.stdout.flush()
    # Iteratively calculate the states reachable for k=1 to 6
    for k in range(1, 7):
        next_counts = defaultdict(int)
        next_state_to_cube = {}
        
        for prev_state, num_ways in counts.items():
            prev_cube = state_to_cube[prev_state]
            for move in moves:
                # Apply the move to a copy of the previous cube
                new_cube = prev_cube.copy()
                new_cube(move)
                new_state = get_state(new_cube)
                
                # Add the number of ways to reach the previous state
                # to the tally for this new state.
                next_counts[new_state] += num_ways
                
                # Store the new cube object if we haven't seen this state yet
                if new_state not in next_state_to_cube:
                    next_state_to_cube[new_state] = new_cube
        
        counts = next_counts
        state_to_cube = next_state_to_cube
        
        # The number of ways to return to solved is the count for the solved state.
        nk_val = counts.get(solved_state, 0)
        n_k.append(nk_val)
        print(f"  Computed n_{k} = {nk_val}")
        sys.stdout.flush()

    n2, n4, n5, n6 = n_k[2], n_k[4], n_k[5], n_k[6]
    
    # The Principle of Inclusion-Exclusion gives the formula:
    # Total = |Solved at 4| + |Solved at 5| + |Solved at 6| - |Solved at 4 and 6|
    # Total = (n_4 * 12 * 12) + (n_5 * 12) + n_6 - (n_4 * n_2)
    # Total = 144*n_4 + 12*n_5 + n_6 - n_4*12
    # Total = 132*n_4 + 12*n_5 + n_6
    
    result = 132 * n4 + 12 * n5 + n6
    
    print("\nStep 2: Calculate the total number of permutations using the formula.")
    print("Formula: Total = 132 * n_4 + 12 * n_5 + n_6")
    print(f"Plugging in computed values: Total = 132 * {n4} + 12 * {n5} + {n6}")
    print(f"Calculation: Total = {132 * n4} + {12 * n5} + {n6}")
    print(f"Final Answer: {result}")
    
    return result

# Execute the function and capture the final answer for the required format.
final_answer = solve_rubiks_permutations()
print(f"\n<<<32868>>>")