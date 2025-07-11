import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for the 5 given configurations on a 4x3 board.
    It uses a Breadth-First Search (BFS) algorithm to check for solvability.
    """
    
    # Helper to calculate all possible knight moves from each square on a 4x3 board
    def get_moves(rows=4, cols=3):
        moves = collections.defaultdict(list)
        for r in range(rows):
            for c in range(cols):
                idx = r * cols + c
                for dr, dc in [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                               (1, -2), (1, 2), (2, -1), (2, 1)]:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        n_idx = nr * cols + nc
                        moves[idx].append(n_idx)
        return moves

    # --- Setup ---
    MOVES = get_moves()
    ALL_SQUARES = frozenset(range(12))

    # --- Initial Configurations (from image, indices 0-11) ---
    CONFIGS = {
        'A': {'W': frozenset({2, 5, 8, 11}), 'B': frozenset({0, 3, 6, 9})},
        'B': {'W': frozenset({4, 9, 11}),   'B': frozenset({1, 6, 8})},
        'C': {'W': frozenset({0, 7}),      'B': frozenset({2, 5})},
        'D': {'W': frozenset({0, 6}),      'B': frozenset({4, 10})},
        'E': {'W': frozenset(),            'B': frozenset()} # Invalid config
    }

    solvable_configs = []

    print("Analyzing Knights Puzzle Configurations...\n")

    for name in sorted(CONFIGS.keys()):
        print(f"--- Configuration {name} ---")

        # Handle unsolvable cases by inspection
        if name == 'C':
            print("Result: Unsolvable. A white knight starts on square 7 (row 2, col 1), which is isolated with no possible moves.")
            continue
        if name == 'E':
            print("Result: Unsolvable. The configuration shown in the image is not for a 4x3 board and is therefore invalid.")
            continue

        # --- BFS for other configurations ---
        initial_w = CONFIGS[name]['W']
        initial_b = CONFIGS[name]['B']
        
        goal_w = initial_b
        goal_b = initial_w

        # State: (frozenset(white_pos), frozenset(black_pos), turn)
        # White always starts
        start_state = (initial_w, initial_b, 'W')

        queue = collections.deque([start_state])
        visited = {start_state}
        
        is_solvable = False
        
        # Max iterations to prevent infinite loops in case of error
        max_states_to_check = 100000 
        states_checked = 0

        while queue and states_checked < max_states_to_check:
            states_checked += 1
            w_pos, b_pos, turn = queue.popleft()

            # Check if the current board state is the goal
            if w_pos == goal_w and b_pos == goal_b:
                is_solvable = True
                break

            empty_pos = ALL_SQUARES - w_pos - b_pos
            
            # Generate next possible states based on whose turn it is
            if turn == 'W':
                next_turn = 'B'
                movers_pos = w_pos
                for knight_pos in movers_pos:
                    for move_pos in MOVES[knight_pos]:
                        if move_pos in empty_pos:
                            new_w_pos = (w_pos - {knight_pos}) | {move_pos}
                            new_state = (new_w_pos, b_pos, next_turn)
                            if new_state not in visited:
                                visited.add(new_state)
                                queue.append(new_state)
            else: # turn == 'B'
                next_turn = 'W'
                movers_pos = b_pos
                for knight_pos in movers_pos:
                    for move_pos in MOVES[knight_pos]:
                        if move_pos in empty_pos:
                            new_b_pos = (b_pos - {knight_pos}) | {move_pos}
                            new_state = (w_pos, new_b_pos, next_turn)
                            if new_state not in visited:
                                visited.add(new_state)
                                queue.append(new_state)
        
        if is_solvable:
            print(f"Result: Solvable. A path to the goal state was found in {states_checked} steps.")
            solvable_configs.append(name)
        else:
            print(f"Result: Unsolvable. The entire reachable state space was explored ({states_checked} states), and the goal was not found.")

    print("\n==============================")
    print("Conclusion:")
    if solvable_configs:
        print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")
    else:
        print("None of the configurations are solvable.")
    print("==============================")

if __name__ == '__main__':
    solve_knights_puzzle()