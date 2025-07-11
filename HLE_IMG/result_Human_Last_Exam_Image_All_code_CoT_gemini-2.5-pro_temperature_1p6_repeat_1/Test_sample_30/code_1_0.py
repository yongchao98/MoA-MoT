import collections

def solve_knights_puzzle():
    """
    Analyzes five configurations of the Knights Puzzle on a 4x3 board to determine which are solvable.
    A puzzle is solvable if the initial positions of white and black knights can be swapped
    through a series of alternating, valid knight moves.
    """

    # The 4x3 board is numbered 0-11 (row-major order).
    # 0  1  2
    # 3  4  5
    # 6  7  8
    # 9 10 11
    # Adjacency list representing all possible knight moves from each square.
    MOVES = {
        0: [5, 7], 1: [6, 8], 2: [3, 7],
        3: [2, 8, 10], 4: [9, 11], 5: [0, 6, 10],
        6: [1, 5, 11], 7: [0, 2], 8: [1, 3],
        9: [4, 8], 10: [3, 5], 11: [4, 6]
    }

    # Initial positions for each configuration, parsed from the provided image.
    # 'w' for white knights, 'b' for black knights.
    configs = {
        'A': {'w': {2, 5, 8, 11}, 'b': {0, 3, 6, 9}},
        'B': {'w': {4, 9, 11},    'b': {1, 6, 8}},
        'C': {'w': {0, 7},        'b': {2, 5}},
        'D': {'w': {0, 7},        'b': {4, 10}},
        'E': {'w': {1, 2, 5},     'b': {0, 3, 4}}
    }

    solvable_configs = []

    for name, pos in configs.items():
        w_initial = tuple(sorted(pos['w']))
        b_initial = tuple(sorted(pos['b']))

        # The final (goal) state has the knight positions swapped.
        w_final = b_initial
        b_final = w_initial

        # A state is defined by (white_positions, black_positions, turn).
        # We start with White's turn.
        initial_state = (w_initial, b_initial, 'W')
        
        # queue for BFS, and a set to store visited states to prevent cycles.
        queue = collections.deque([initial_state])
        visited = {initial_state}
        
        is_solvable = False
        
        # Use a limit to prevent extremely long searches for complex unsolvable cases.
        max_states_to_check = 200000 
        states_checked = 0

        while queue:
            states_checked += 1
            if states_checked > max_states_to_check:
                break # Assume unsolvable if search is too long.

            w_curr, b_curr, turn = queue.popleft()

            # Check if the current state is the goal configuration.
            if w_curr == w_final and b_curr == b_final:
                is_solvable = True
                break

            occupied_squares = set(w_curr) | set(b_curr)

            if turn == 'W':
                # Generate next states from White's possible moves.
                for i, p_w in enumerate(w_curr):
                    for move in MOVES.get(p_w, []):
                        if move not in occupied_squares:
                            w_next_list = list(w_curr)
                            w_next_list[i] = move
                            w_next = tuple(sorted(w_next_list))
                            
                            next_state = (w_next, b_curr, 'B')
                            if next_state not in visited:
                                visited.add(next_state)
                                queue.append(next_state)
            else: # turn == 'B'
                # Generate next states from Black's possible moves.
                for i, p_b in enumerate(b_curr):
                    for move in MOVES.get(p_b, []):
                        if move not in occupied_squares:
                            b_next_list = list(b_curr)
                            b_next_list[i] = move
                            b_next = tuple(sorted(b_next_list))
                            
                            next_state = (w_curr, b_next, 'W')
                            if next_state not in visited:
                                visited.add(next_state)
                                queue.append(next_state)
        
        if is_solvable:
            solvable_configs.append(name)
            print(f"Configuration {name}: Solvable")
        else:
            print(f"Configuration {name}: Unsolvable")

    print("\n--- Summary ---")
    if solvable_configs:
        print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")
    else:
        print("None of the configurations are solvable.")

# Run the solver.
solve_knights_puzzle()