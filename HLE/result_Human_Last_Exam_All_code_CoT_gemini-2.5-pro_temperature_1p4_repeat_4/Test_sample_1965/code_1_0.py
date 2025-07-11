# This script requires the 'pycuber' library.
# You can install it by running: pip install pycuber

import pycuber as pc
from collections import defaultdict

def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at the 4th, 5th, or 6th turn.
    """

    # The 12 standard 90-degree moves
    moves = [
        pc.Move(m) for m in ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]
    ]

    # Get the tuple representation of a solved cube state, which is hashable.
    solved_state = pc.Cube().get_cube()

    # W_k maps a state (tuple) to the number of ways to reach it in k moves.
    W_k = {solved_state: 1}
    
    # N[k] will store the number of k-move sequences that result in the solved state.
    N = {0: 1}

    max_k = 6
    for k in range(1, max_k + 1):
        W_k_plus_1 = defaultdict(int)
        
        # For each state reachable in k moves...
        for state_tuple, count in W_k.items():
            # Create a cube object representing this state to apply new moves.
            temp_cube = pc.Cube(state_tuple)
            
            # Apply each of the 12 possible moves.
            for move in moves:
                # Create a copy to avoid modifying the cube in the loop.
                moved_cube = temp_cube.copy()
                moved_cube.do(move)
                new_state_tuple = moved_cube.get_cube()
                
                # Add the number of paths to the new state.
                W_k_plus_1[new_state_tuple] += count
        
        W_k = W_k_plus_1
        # Store the number of ways to return to the solved state.
        N[k] = W_k.get(solved_state, 0)

    # Values needed for the formula
    n2 = N.get(2, 0)
    n4 = N.get(4, 0)
    n5 = N.get(5, 0)
    n6 = N.get(6, 0)
    
    # Using the Inclusion-Exclusion principle:
    # Total = |A| + |B| + |C| - |A ∩ C|
    # Total = (N(4) * 12^2) + (N(5) * 12) + N(6) - (N(4) * N(2))
    
    term_A = n4 * 144
    term_B = n5 * 12
    term_C = n6
    term_AC = n4 * n2

    total_permutations = term_A + term_B + term_C - term_AC

    print("Calculation based on computed N(k) values:")
    print(f"N(2): Number of 2-move sequences returning to solved = {n2}")
    print(f"N(4): Number of 4-move sequences returning to solved = {n4}")
    print(f"N(5): Number of 5-move sequences returning to solved = {n5}")
    print(f"N(6): Number of 6-move sequences returning to solved = {n6}")
    print("\nApplying the Principle of Inclusion-Exclusion:")
    print(f"Total = (N(4) * 144) + (N(5) * 12) + N(6) - (N(4) * N(2))")
    print(f"Total = ({n4} * 144) + ({n5} * 12) + {n6} - ({n4} * {n2})")
    print(f"Total = {term_A} + {term_B} + {term_C} - {term_AC}")
    print(f"\nFinal Answer: {total_permutations}")

solve_rubiks_permutations()
<<<19056>>>