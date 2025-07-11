import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for five configurations on a 4x3 board.

    The state of the puzzle is represented by a tuple of length 12 (for the
    12 squares) and the current player's turn ('W' for white, 'B' for black).
    A Breadth-First Search (BFS) is used to explore the state space from
    each initial configuration. A configuration is solvable if the goal state
    (where white and black knights have swapped positions) is reachable.
    """

    # The 4x3 board is numbered 0 to 11, row by row:
    #  0  1  2
    #  3  4  5
    #  6  7  8
    #  9 10 11
    # Pre-calculated valid knight moves for each square.
    MOVES = {
        0: [5, 7], 1: [6, 8], 2: [3, 7],
        3: [2, 8, 10], 4: [9, 11], 5: [0, 6, 10],
        6: [1, 5, 11], 7: [0, 2], 8: [1, 3],
        9: [4, 8], 10: [3, 5], 11: [4, 6]
    }

    # Initial positions of black and white knights for configurations A-E.
    # The keys 'initial_B' and 'initial_W' hold sets of positions.
    configs = {
        'A': {'initial_B': {0, 3, 6, 9}, 'initial_W': {2, 5, 8, 11}},
        'B': {'initial_B': {1, 6, 8},   'initial_W': {3, 10, 11}},
        'C': {'initial_B': {2, 5},      'initial_W': {0, 7}},
        'D': {'initial_B': {4, 9},      'initial_W': {1, 8}},
        'E': {'initial_B': {0, 3, 4},   'initial_W': {1, 2, 5}}
    }

    def create_board_tuple(black_knights, white_knights, size=12):
        """Helper to create a board tuple from knight positions."""
        board = ['E'] * size
        for pos in black_knights:
            board[pos] = 'B'
        for pos in white_knights:
            board[pos] = 'W'
        return tuple(board)

    solvable_configs = []

    print("Analyzing configurations...")
    # Iterate through each configuration to check for solvability.
    for name, config_data in sorted(configs.items()):
        initial_w_pos = config_data['initial_W']
        initial_b_pos = config_data['initial_B']

        # The initial board state.
        initial_board = create_board_tuple(initial_b_pos, initial_w_pos)
        # The goal state is where knights of opposite colors have swapped places.
        goal_board = create_board_tuple(initial_w_pos, initial_b_pos)

        # Initialize BFS queue and visited set. The state includes the turn.
        # White always moves first.
        start_state = (initial_board, 'W')
        queue = collections.deque([start_state])
        visited = {start_state}
        
        is_solvable = False
        
        # A safety limit on iterations to prevent extremely long runs.
        max_iterations = 250000
        count = 0
        
        while queue:
            if count > max_iterations:
                print(f"Configuration {name} search limit exceeded.")
                break
            count += 1
            
            current_board, current_turn = queue.popleft()

            if current_board == goal_board:
                is_solvable = True
                break

            piece_to_move = current_turn
            next_turn = 'B' if current_turn == 'W' else 'W'
            
            # Find all possible moves for the current player.
            for from_pos in range(12):
                if current_board[from_pos] == piece_to_move:
                    # Check all valid knight moves from this square.
                    for to_pos in MOVES[from_pos]:
                        # A knight can only move to an empty square.
                        if current_board[to_pos] == 'E':
                            # Generate the board for the next state.
                            board_list = list(current_board)
                            board_list[from_pos], board_list[to_pos] = board_list[to_pos], board_list[from_pos]
                            next_board = tuple(board_list)
                            
                            next_state = (next_board, next_turn)
                            
                            # If we haven't seen this state before, add it to the queue.
                            if next_state not in visited:
                                visited.add(next_state)
                                queue.append(next_state)
        
        if is_solvable:
            print(f"Configuration {name} is solvable.")
            solvable_configs.append(name)
        else:
            print(f"Configuration {name} is not solvable.")

    # Print the final result.
    final_answer = ", ".join(sorted(solvable_configs))
    print("\nConclusion:")
    print(f"The solvable configurations are: {final_answer}")

# Execute the solver function
solve_knights_puzzle()