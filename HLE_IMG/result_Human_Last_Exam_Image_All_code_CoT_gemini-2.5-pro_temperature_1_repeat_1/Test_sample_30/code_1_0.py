import collections

def solve_knights_puzzle():
    """
    This script determines which of the five knight puzzle configurations are solvable.
    It models the puzzle as a state-space search and uses Breadth-First Search (BFS)
    to find if a solution exists.
    """
    
    # Board dimensions
    ROWS, COLS = 4, 3

    # Pre-calculate valid knight moves for each square to speed up the search
    # (r, c) coordinates are used for board positions.
    moves_map = {}
    for r in range(ROWS):
        for c in range(COLS):
            moves_map[(r, c)] = []
            # All 8 possible L-shaped moves for a knight
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                           (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < ROWS and 0 <= nc < COLS:
                    moves_map[(r, c)].append((nr, nc))

    def solve_config(initial_board):
        """
        Solves the knight puzzle for a given initial configuration using BFS.

        Args:
            initial_board: A tuple of tuples representing the initial board state.
                           1 for white, -1 for black, 0 for empty.

        Returns:
            True if solvable, False otherwise.
        """
        # The goal board is created by swapping the colors of all pieces
        goal_board = tuple(tuple(-piece for piece in row) for row in initial_board)

        # The state is a tuple: (board_configuration, player_turn)
        # Player turn: 1 for White, -1 for Black. White starts.
        start_state = (initial_board, 1)

        # Queue for BFS, initialized with the starting state
        queue = collections.deque([start_state])
        # A set to keep track of visited states to avoid cycles and redundant work
        visited = {start_state}

        while queue:
            current_board, turn = queue.popleft()

            # If the current board matches the goal, a solution is found
            if current_board == goal_board:
                return True

            # Find all pieces for the current player
            piece_positions = []
            for r in range(ROWS):
                for c in range(COLS):
                    if current_board[r][c] == turn:
                        piece_positions.append((r, c))

            # Generate next possible states from the current state
            for r_from, c_from in piece_positions:
                for r_to, c_to in moves_map[(r_from, c_from)]:
                    # A move is valid only if the destination square is empty
                    if current_board[r_to][c_to] == 0:
                        # Create the new board by moving the piece
                        new_board_list = [list(row) for row in current_board]
                        new_board_list[r_to][c_to] = turn
                        new_board_list[r_from][c_from] = 0
                        new_board = tuple(tuple(row) for row in new_board_list)
                        
                        # Create the new state with the turn switched
                        next_turn = -turn
                        new_state = (new_board, next_turn)

                        if new_state not in visited:
                            visited.add(new_state)
                            queue.append(new_state)

        # If the queue becomes empty and the goal was not reached, it's unsolvable
        return False

    # Define the 5 initial configurations based on the image
    # 1: White Knight, -1: Black Knight, 0: Empty Square
    configurations = {
        'A': ((-1, 0, 1),
              (-1, 0, 1),
              (-1, 0, 1),
              (-1, 0, 1)),
        'B': (( 0,-1, 0),
              ( 0, 1, 0),
              (-1, 0,-1),
              ( 0, 1, 1)),
        'C': (( 1, 0,-1),
              ( 0, 0,-1),
              ( 0, 1, 0),
              ( 0, 0, 0)),
        'D': (( 0, 1, 0),
              (-1, 0, 0),
              ( 0, 1, 0),
              (-1, 0, 0)),
        'E': ((-1, 1, 1),
              (-1,-1, 1),
              ( 0, 0, 0),
              ( 0, 0, 0))
    }
    
    solvable_configs = []
    print("Analyzing the Knight Puzzle configurations...")
    for name, board in configurations.items():
        is_solvable = solve_config(board)
        if is_solvable:
            solvable_configs.append(name)
    
    print("\nResult:")
    if solvable_configs:
        print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")
    else:
        print("None of the provided configurations are solvable.")

solve_knights_puzzle()