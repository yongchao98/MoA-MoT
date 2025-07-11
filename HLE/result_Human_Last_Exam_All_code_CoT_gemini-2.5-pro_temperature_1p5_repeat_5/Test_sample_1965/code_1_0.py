import sys
from collections import defaultdict

def solve_rubiks_permutations():
    """
    This script calculates the number of 6-move sequences that solve a Rubik's cube
    at some point during the final 3 moves.
    """

    # --- Cube representation and move definitions ---
    
    # A cube state is a tuple of 54 sticker values.
    # Solved state: 0s for U face, 1s for D face, etc.
    SOLVED_STATE = tuple(c for c in "UUUUUUUUUDDDDDDDDDLLLLLLLLLRRRRRRRRRFFFFFFFFFBBBBBBBBB")
    
    # Pre-computed permutation maps for each of the 6 clockwise 90-degree moves.
    # Each map is a list of 54 indices, where map[i] is the source index for the sticker at new position i.
    # E.g., new_state[i] = old_state[map[i]]
    MOVE_PERMS = {
        'U': (2, 5, 8, 1, 4, 7, 0, 3, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
              36, 37, 38, 21, 22, 23, 24, 25, 26, 18, 19, 20, 30, 31, 32, 33, 34, 35, 
              45, 46, 47, 39, 40, 41, 42, 43, 44, 27, 28, 29, 48, 49, 50, 51, 52, 53),
        'D': (0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 14, 17, 10, 13, 16, 9, 12, 15,
              18, 19, 20, 21, 22, 23, 44, 43, 42, 27, 28, 29, 30, 31, 32, 53, 52, 51,
              36, 37, 38, 39, 40, 41, 26, 25, 24, 45, 46, 47, 48, 49, 50, 34, 33, 32),
        'L': (0, 1, 47, 3, 4, 50, 6, 7, 53, 42, 10, 11, 39, 13, 14, 36, 16, 17,
              20, 23, 26, 19, 22, 25, 18, 21, 24, 27, 28, 2, 30, 31, 5, 33, 34, 8,
              15, 37, 38, 12, 40, 41, 9, 43, 44, 45, 46, 49, 48, 1, 51, 52, 3),
        'R': (27, 1, 2, 30, 4, 5, 33, 7, 8, 9, 10, 44, 12, 13, 41, 15, 16, 38,
              18, 19, 20, 21, 22, 23, 24, 25, 26, 29, 32, 35, 28, 31, 34, 27, 30, 33,
              36, 37, 17, 39, 40, 14, 42, 43, 11, 6, 46, 47, 4, 49, 50, 2, 52, 53, 0),
        'F': (0, 1, 2, 3, 4, 5, 20, 23, 26, 11, 10, 9, 12, 13, 14, 15, 16, 17,
              18, 19, 33, 21, 22, 30, 24, 25, 27, 6, 28, 29, 7, 31, 32, 8, 34, 35,
              38, 41, 44, 37, 40, 43, 36, 39, 42, 45, 46, 47, 48, 49, 50, 51, 52, 53),
        'B': (51, 52, 53, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 29, 32, 35,
              2, 19, 20, 5, 22, 23, 4, 25, 26, 27, 28, 15, 30, 31, 16, 33, 34, 17,
              36, 37, 38, 39, 40, 41, 42, 43, 44, 47, 50, 45, 46, 49, 48, 18, 21, 24)
    }

    def invert_perm(p):
        inv = [0] * len(p)
        for i, val in enumerate(p):
            inv[val] = i
        return tuple(inv)

    # Generate all 12 move permutations (6 clockwise, 6 anti-clockwise)
    all_moves = []
    for name, p_map in MOVE_PERMS.items():
        all_moves.append(p_map)  # Clockwise move
        all_moves.append(invert_perm(p_map)) # Anti-clockwise move

    def apply_move(state, perm):
        return tuple(state[i] for i in perm)

    # --- Path counting algorithm ---
    
    # counts[k] stores {state: num_paths} for paths of length k
    counts = {SOLVED_STATE: 1}
    n = {} # Stores n(k) = number of ways to return to SOLVED_STATE in k moves
    max_k = 6
    
    for k in range(1, max_k + 1):
        next_counts = defaultdict(int)
        for state, num_paths in counts.items():
            for move_perm in all_moves:
                next_state = apply_move(state, move_perm)
                next_counts[next_state] += num_paths
        counts = next_counts
        n[k] = counts.get(SOLVED_STATE, 0)

    # --- Calculate and print the final answer ---
    
    # As derived, Total = 132 * n(4) + n(6)
    n4 = n.get(4, 0)
    n6 = n.get(6, 0)
    
    result = 132 * n4 + n6
    
    print(f"The number of ways to return to the solved state in exactly 4 moves is n(4) = {n4}.")
    print(f"The number of ways to return to the solved state in exactly 6 moves is n(6) = {n6}.")
    print("\nThe total number of permutations that result in the cube returning to its original configuration at some point during the final 3 moves is given by the formula 132 * n(4) + n(6).")
    print(f"\nCalculation: 132 * {n4} + {n6} = {result}")

solve_rubiks_permutations()