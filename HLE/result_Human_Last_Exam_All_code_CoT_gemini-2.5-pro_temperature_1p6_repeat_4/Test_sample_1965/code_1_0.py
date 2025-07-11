import pycuber as pc
from collections import defaultdict
import sys

def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at step 4, 5, or 6.
    """
    # A hashable representation of a cube's state.
    def get_state_key(cube):
        return frozenset(c for c in cube.cubies)

    # The 12 standard 90-degree moves.
    moves = [pc.Move(s) for s in ("U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'")]

    # We start with one way to have a solved cube in 0 moves (the empty sequence).
    # reachable_states maps a state's key to its cube object and the number of ways to reach it.
    solved_cube = pc.Cube()
    solved_key = get_state_key(solved_cube)
    reachable_states = {solved_key: (solved_cube, 1)}

    k_is_I = {}

    print("Thinking Process:")
    print("Finding the number of move sequences of length n that result in a solved cube...")
    # Explore the state space step by step for 6 moves.
    for n in range(1, 7):
        next_reachable_states = defaultdict(lambda: [None, 0])
        
        for cube_obj, count in reachable_states.values():
            for move in moves:
                next_cube = cube_obj.copy()
                next_cube.perform_move(move)
                next_key = get_state_key(next_cube)
                
                # Aggregate counts for states reached via different paths.
                # The cube object stored is arbitrary; we just need one for the next expansion.
                entry = next_reachable_states[next_key]
                entry[0] = next_cube
                entry[1] += count
        
        # Convert defaultdict back to dict for the next iteration.
        reachable_states = {k: (v[0], v[1]) for k, v in next_reachable_states.items()}
        
        # Store the number of ways to return to the solved state for step n.
        k_is_I[n] = reachable_states.get(solved_key, (None, 0))[1]
        print(f"k_{n}_is_I (Ways to be solved in {n} moves): {k_is_I[n]}")

    # Extract the required values.
    k_4_is_I = k_is_I[4]
    k_5_is_I = k_is_I[5]
    k_6_is_I = k_is_I[6]
    
    # Calculate the total using the derived formula.
    total = 132 * k_4_is_I + 12 * k_5_is_I + k_6_is_I
    
    print("\nFinal Calculation:")
    print(f"The number of desired permutations is given by the formula:")
    print(f"Total = 132 * k_4_is_I + 12 * k_5_is_I + k_6_is_I")
    print(f"Substituting the computed values:")
    print(f"Total = 132 * {k_4_is_I} + 12 * {k_5_is_I} + {k_6_is_I}")
    print(f"Total = {132 * k_4_is_I} + {12 * k_5_is_I} + {k_6_is_I}")
    print(f"Total = {total}")

solve_rubiks_permutations()
<<<42912>>>