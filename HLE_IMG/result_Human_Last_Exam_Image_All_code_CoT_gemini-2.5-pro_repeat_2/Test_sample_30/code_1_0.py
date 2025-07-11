import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for five given initial configurations on a 4x3 board.
    It determines which configurations are solvable, i.e., where white and black
    knights can swap positions through a series of alternating moves.
    """
    
    # --- 1. Setup Board and Moves ---
    BOARD_ROWS = 4
    BOARD_COLS = 3
    BOARD_SIZE = BOARD_ROWS * BOARD_COLS

    def get_all_knight_moves():
        """Pre-computes all possible knight moves for each square on the 4x3 board."""
        all_moves = collections.defaultdict(list)
        for r in range(BOARD_ROWS):
            for c in range(BOARD_COLS):
                pos = r * BOARD_COLS + c
                deltas = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                          (1, -2), (1, 2), (2, -1), (2, 1)]
                for dr, dc in deltas:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < BOARD_ROWS and 0 <= nc < BOARD_COLS:
                        all_moves[pos].append(nr * BOARD_COLS + nc)
        return all_moves

    ALL_MOVES = get_all_knight_moves()

    def board_from_positions(white_pos, black_pos):
        """Creates a board tuple from lists of knight positions."""
        board = ['E'] * BOARD_SIZE
        for pos in white_pos:
            board[pos] = 'W'
        for pos in black_pos:
            board[pos] = 'B'
        return tuple(board)

    # --- 2. BFS Solver ---
    def is_solvable(initial_board, goal_board):
        """
        Determines if a puzzle is solvable using Breadth-First Search (BFS).
        A state is a tuple: (board_configuration, player_to_move).
        """
        queue = collections.deque([(initial_board, 'W')])
        visited = {(initial_board, 'W')}

        while queue:
            current_board, current_player = queue.popleft()

            if current_board == goal_board:
                return True

            piece_to_move = current_player
            next_player = 'B' if current_player == 'W' else 'W'
            
            for i, piece in enumerate(current_board):
                if piece == piece_to_move:
                    start_pos = i
                    for end_pos in ALL_MOVES[start_pos]:
                        if current_board[end_pos] == 'E':
                            new_board_list = list(current_board)
                            new_board_list[start_pos], new_board_list[end_pos] = 'E', piece_to_move
                            new_board = tuple(new_board_list)
                            
                            new_state = (new_board, next_player)
                            if new_state not in visited:
                                visited.add(new_state)
                                queue.append(new_state)
        return False

    # --- 3. Define Configurations and Run Solver ---
    # Positions are 0-indexed: 0-2 (row 1), 3-5 (row 2), etc.
    # pos = row * 3 + col
    configs = {
        'A': {'W': [2, 5, 8, 11], 'B': [0, 3, 6, 9]},
        'B': {'W': [4, 9, 11], 'B': [1, 6, 8]},
        'C': {'W': [0, 7], 'B': [2, 5]},
        'D': {'W': [0, 7], 'B': [4, 10]},
        'E': {'W': [1, 2, 5], 'B': [0, 3, 4]}
    }
    
    solvable_configs = []
    
    print("Analyzing Knights Puzzle configurations...")
    for name, positions in configs.items():
        initial_board = board_from_positions(positions['W'], positions['B'])
        goal_board = board_from_positions(positions['B'], positions['W'])
        
        if is_solvable(initial_board, goal_board):
            solvable_configs.append(name)
    
    print("\n--- Results ---")
    if solvable_configs:
        print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")
    else:
        print("None of the configurations were found to be solvable.")

solve_knights_puzzle()