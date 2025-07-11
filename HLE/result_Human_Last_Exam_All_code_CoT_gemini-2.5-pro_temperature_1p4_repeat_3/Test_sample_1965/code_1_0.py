import collections

def get_moves():
    """
    Defines the permutations for the 12 standard 90-degree moves.
    Each move is a list of cycles of sticker indices.
    U:0-8, L:9-17, F:18-26, R:27-35, B:36-44, D:45-53
    """
    # Clockwise 90-degree turns
    moves_cw = {
        'U': [(0, 2, 8, 6), (1, 5, 7, 3), (9, 18, 27, 36), (10, 19, 28, 37), (11, 20, 29, 38)],
        'D': [(45, 47, 53, 51), (46, 50, 52, 48), (24, 15, 42, 33), (25, 16, 43, 34), (26, 17, 44, 35)],
        'L': [(9, 11, 17, 15), (10, 14, 16, 12), (0, 44, 45, 18), (3, 41, 48, 21), (6, 38, 51, 24)],
        'R': [(27, 29, 35, 33), (28, 32, 34, 30), (2, 20, 47, 42), (5, 23, 50, 39), (8, 26, 53, 36)],
        'F': [(18, 20, 26, 24), (19, 23, 25, 21), (6, 11, 47, 33), (7, 14, 46, 30), (8, 17, 45, 27)],
        'B': [(36, 38, 44, 42), (37, 41, 43, 39), (0, 29, 53, 15), (1, 32, 52, 12), (2, 35, 51, 9)]
    }

    # Generate counter-clockwise moves by reversing the cycles
    moves_ccw = { key + "'": [cycle[::-1] for cycle in val] for key, val in moves_cw.items() }
    
    # Combine all 12 moves
    all_moves = list(moves_cw.values()) + list(moves_ccw.values())
    return all_moves

def apply_permutation(state_tuple, permutation):
    """Applies a permutation to a given state tuple."""
    state_list = list(state_tuple)
    new_state_list = list(state_tuple)
    for cycle in permutation:
        for i in range(len(cycle)):
            new_state_list[cycle[i]] = state_list[cycle[i-1]]
    return tuple(new_state_list)

def solve_rubiks_problem():
    """
    Calculates the number of permutations that solve the cube during the final 3 moves.
    """
    print("Starting calculation... This may take a few seconds.")

    # Step 1: Get move definitions and set up initial state
    moves = get_moves()
    # Initial state: 9 stickers of color 0, then 9 of color 1, etc.
    initial_state = tuple(i // 9 for i in range(54))
    
    # counts will store {state: number_of_ways_to_reach_it}
    counts = {initial_state: 1}
    n_values = {}

    # Step 2: Use BFS to find N(k) for k=1 to 6
    for k in range(1, 7):
        next_counts = collections.defaultdict(int)
        for state, num_ways in counts.items():
            for move in moves:
                next_state = apply_permutation(state, move)
                next_counts[next_state] += num_ways
        counts = next_counts
        # N(k) is the number of ways to return to the initial state
        n_values[k] = counts.get(initial_state, 0)
        # print(f"Computed N({k}) = {n_values[k]}") # For debugging

    N4 = n_values[4]
    N5 = n_values[5]
    N6 = n_values[6]

    # Step 3: Apply the inclusion-exclusion formula
    result = 132 * N4 + 12 * N5 + N6

    # Step 4: Print the detailed calculation and the final answer
    print("\nBased on the Principle of Inclusion-Exclusion, the formula is: 132 * N(4) + 12 * N(5) + N(6)")
    print("\nWhere N(k) is the number of ways to return to the solved state in exactly k moves.")
    print(f"Number of ways for 4 moves (N(4)): {N4}")
    print(f"Number of ways for 5 moves (N(5)): {N5}")
    print(f"Number of ways for 6 moves (N(6)): {N6}")
    
    print("\nFinal Calculation:")
    print(f"Result = 132 * {N4} + 12 * {N5} + {N6} = {result}")

solve_rubiks_problem()
<<<18936>>>