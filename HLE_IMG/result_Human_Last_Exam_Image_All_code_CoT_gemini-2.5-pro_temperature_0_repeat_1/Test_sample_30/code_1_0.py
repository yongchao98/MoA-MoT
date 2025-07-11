import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for five given configurations on a 4x3 board.

    The method uses Breadth-First Search (BFS) to explore the state space
    of the puzzle. A state is defined by the arrangement of knights on the
    board and which color's turn it is to move. The goal is to determine
    if the initial positions of white and black knights can be swapped.
    """

    # Board dimensions and setup
    WIDTH = 3
    HEIGHT = 4
    NUM_SQUARES = WIDTH * HEIGHT

    def get_knight_moves():
        """Pre-calculates all possible knight moves from each square."""
        moves = collections.defaultdict(list)
        for r in range(HEIGHT):
            for c in range(WIDTH):
                idx = r * WIDTH + c
                for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                               (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < HEIGHT and 0 <= nc < WIDTH:
                        n_idx = nr * WIDTH + nc
                        moves[idx].append(n_idx)
        return moves

    KNIGHT_MOVES = get_knight_moves()

    def is_solvable(initial_config):
        """
        Determines if a configuration is solvable using BFS.
        """
        # 1. Define initial and goal board layouts.
        # A board is a tuple of 12 characters ('W', 'B', 'E').
        initial_board_list = ['E'] * NUM_SQUARES
        for pos in initial_config['W']:
            initial_board_list[pos] = 'W'
        for pos in initial_config['B']:
            initial_board_list[pos] = 'B'
        initial_board = tuple(initial_board_list)

        goal_board_list = ['E'] * NUM_SQUARES
        for pos in initial_config['B']:
            goal_board_list[pos] = 'W'
        for pos in initial_config['W']:
            goal_board_list[pos] = 'B'
        goal_board = tuple(goal_board_list)

        # 2. Initialize BFS queue and visited set.
        # A state is (board_tuple, next_player_to_move). White starts.
        start_state = (initial_board, 'W')
        queue = collections.deque([start_state])
        visited = {start_state}

        # 3. Perform BFS.
        while queue:
            current_board, current_player = queue.popleft()

            # Check if the goal configuration is reached.
            if current_board == goal_board:
                return True

            # 4. Generate next states.
            next_player = 'B' if current_player == 'W' else 'W'
            for i in range(NUM_SQUARES):
                if current_board[i] == current_player:
                    # Try all moves for the current knight.
                    for move_to_idx in KNIGHT_MOVES[i]:
                        if current_board[move_to_idx] == 'E':
                            # Create the new board state.
                            new_board_list = list(current_board)
                            new_board_list[i], new_board_list[move_to_idx] = 'E', current_player
                            new_board = tuple(new_board_list)
                            
                            new_state = (new_board, next_player)
                            if new_state not in visited:
                                visited.add(new_state)
                                queue.append(new_state)
        
        # If the queue is empty and goal not found, it's unsolvable.
        return False

    # Define the 5 initial configurations from the image.
    # Board indexed 0-11, row by row.
    # 0  1  2
    # 3  4  5
    # 6  7  8
    # 9 10 11
    configs = {
        'A': {'B': [0, 3, 6, 9], 'W': [2, 5, 8, 11]},
        'B': {'B': [0, 3, 6, 8], 'W': [1, 5, 9, 11]},
        'C': {'B': [2, 5], 'W': [0, 8]},
        'D': {'B': [2, 5, 11], 'W': [0, 6, 8]},
        'E': {'B': [0, 3, 5], 'W': [1, 2, 8]}
    }

    print("Analyzing the Knights Puzzle configurations...")
    solvable_configs = []
    # Sort keys to ensure consistent output order
    for name in sorted(configs.keys()):
        config = configs[name]
        if is_solvable(config):
            solvable_configs.append(name)
            print(f"Result: Configuration {name} is SOLVABLE.")
        else:
            print(f"Result: Configuration {name} is UNSOLVABLE.")
    
    print("\n--- Conclusion ---")
    if not solvable_configs:
        print("None of the five configurations are solvable.")
    else:
        print(f"The solvable configuration(s) are: {', '.join(solvable_configs)}")

# Run the solver
solve_knights_puzzle()