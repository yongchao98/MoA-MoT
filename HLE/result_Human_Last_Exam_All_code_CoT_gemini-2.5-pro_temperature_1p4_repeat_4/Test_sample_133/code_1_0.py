import copy

def solve_connect_four():
    """
    Analyzes the given Connect 4 board state to find the optimal moves for player 'O'.
    An optimal move is one that leads to a win in the minimum number of turns.
    """

    # --- Setup ---
    ROWS = 6
    COLS = 7
    PLAYER_O = 'O'
    PLAYER_X = 'X'
    EMPTY = '.'

    # The board state as given in the problem
    board = [
        [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY],
        [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY],
        [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY],
        [EMPTY, EMPTY, EMPTY, PLAYER_O, PLAYER_O, EMPTY, EMPTY],
        [PLAYER_O, EMPTY, PLAYER_X, PLAYER_O, PLAYER_X, PLAYER_X, PLAYER_X],
        [PLAYER_X, PLAYER_O, PLAYER_O, PLAYER_X, PLAYER_X, PLAYER_O, PLAYER_X],
    ]

    # --- Helper Functions ---
    def get_valid_moves(b):
        """Returns a list of columns that are not full."""
        return [c for c in range(COLS) if b[0][c] == EMPTY]

    def get_next_open_row(b, col):
        """Returns the row index where a piece would fall."""
        for r in range(ROWS - 1, -1, -1):
            if b[r][col] == EMPTY:
                return r
        return None

    def drop_piece(b, col, piece):
        """Returns a new board state after dropping a piece."""
        b_copy = copy.deepcopy(b)
        row = get_next_open_row(b_copy, col)
        if row is not None:
            b_copy[row][col] = piece
        return b_copy

    def check_win(b, piece):
        """Checks if the specified player has won."""
        # Horizontal check
        for r in range(ROWS):
            for c in range(COLS - 3):
                if all(b[r][c+i] == piece for i in range(4)):
                    return True
        # Vertical check
        for c in range(COLS):
            for r in range(ROWS - 3):
                if all(b[r+i][c] == piece for i in range(4)):
                    return True
        # Diagonal (down-right) check
        for r in range(ROWS - 3):
            for c in range(COLS - 3):
                if all(b[r+i][c+i] == piece for i in range(4)):
                    return True
        # Diagonal (up-right) check
        for r in range(3, ROWS):
            for c in range(COLS - 3):
                if all(b[r-i][c+i] == piece for i in range(4)):
                    return True
        return False

    # --- Main Logic ---

    optimal_moves = []
    
    # First, check for any immediate winning moves (win in 1)
    for move_col in get_valid_moves(board):
        temp_board = drop_piece(board, move_col, PLAYER_O)
        if check_win(temp_board, PLAYER_O):
            optimal_moves.append(move_col)
    
    # If no win in 1, search for a win in 2 (a move creating a fork)
    if not optimal_moves:
        for move_col in get_valid_moves(board):
            board_after_o = drop_piece(board, move_col, PLAYER_O)
            
            # This move is a "fork" if for ALL opponent replies, 'O' can win.
            is_a_fork = True
            opponent_moves = get_valid_moves(board_after_o)
            
            # If the opponent has no moves, it's a draw, not a win.
            if not opponent_moves:
                is_a_fork = False

            for opp_col in opponent_moves:
                board_after_x = drop_piece(board_after_o, opp_col, PLAYER_X)
                
                # Check if 'O' has at least one winning move now
                has_o_win = False
                for o_final_col in get_valid_moves(board_after_x):
                    if check_win(drop_piece(board_after_x, o_final_col, PLAYER_O), PLAYER_O):
                        has_o_win = True
                        break
                
                if not has_o_win:
                    # Opponent found a reply to stop the win, so it wasn't a fork.
                    is_a_fork = False
                    break
            
            if is_a_fork:
                optimal_moves.append(move_col)

    # --- Format and Print Output ---
    if not optimal_moves:
        print("No optimal winning moves found.")
    else:
        final_moves_str = []
        col_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        
        optimal_moves.sort()
        
        for col_idx in optimal_moves:
            row_idx = get_next_open_row(board, col_idx)
            row_notation = ROWS - row_idx
            final_moves_str.append(f"{col_names[col_idx]}{row_notation}")
            
        print(", ".join(final_moves_str))

solve_connect_four()