import copy

def solve_connect4():
    """
    Finds the optimal moves for 'O' in the given Connect 4 game state.
    """
    # Board representation: 6 rows, 7 columns
    # 0,0 is top-left
    # Rows 0-5 correspond to board rows 1-6
    # Cols 0-6 correspond to board columns a-g
    board = [
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', 'O', '.', '.'],
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    ]
    rows = 6
    cols = 7
    player = 'O'
    opponent = 'X'

    def get_valid_moves(current_board):
        """Returns a list of columns that are not full."""
        return [c for c in range(cols) if current_board[0][c] == '.']

    def get_next_open_row(current_board, col):
        """Returns the row an object would fall to in a given column."""
        for r in range(rows - 1, -1, -1):
            if current_board[r][col] == '.':
                return r
        return -1

    def check_win(current_board, p):
        """Checks if player p has won."""
        # Check horizontal
        for r in range(rows):
            for c in range(cols - 3):
                if all(current_board[r][c + i] == p for i in range(4)):
                    return True
        # Check vertical
        for c in range(cols):
            for r in range(rows - 3):
                if all(current_board[r + i][c] == p for i in range(4)):
                    return True
        # Check diagonal \
        for r in range(rows - 3):
            for c in range(cols - 3):
                if all(current_board[r + i][c + i] == p for i in range(4)):
                    return True
        # Check diagonal /
        for r in range(3, rows):
            for c in range(cols - 3):
                if all(current_board[r - i][c + i] == p for i in range(4)):
                    return True
        return False

    def col_to_char(c):
        return chr(ord('a') + c)
        
    def row_to_char(r):
        return str(r + 1)

    # --- Main Logic ---

    # Level 1: Find any win in 1 move.
    win_in_1_moves = []
    for move_col in get_valid_moves(board):
        temp_board = copy.deepcopy(board)
        row = get_next_open_row(temp_board, move_col)
        temp_board[row][move_col] = player
        if check_win(temp_board, player):
            win_in_1_moves.append(f"{col_to_char(move_col)}{row_to_char(row)}")
    
    if win_in_1_moves:
        print(', '.join(sorted(win_in_1_moves)))
        return

    # Level 2: Find any win in 2 moves (forced win).
    win_in_2_moves = []
    for move_col in get_valid_moves(board):
        board_after_O1 = copy.deepcopy(board)
        row1 = get_next_open_row(board_after_O1, move_col)
        board_after_O1[row1][move_col] = player

        is_forced_win = True
        opponent_moves = get_valid_moves(board_after_O1)
        
        # If board is full after O's move, it's a draw, not a win.
        if not opponent_moves:
            is_forced_win = False
            
        for opp_move_col in opponent_moves:
            board_after_X = copy.deepcopy(board_after_O1)
            row_x = get_next_open_row(board_after_X, opp_move_col)
            board_after_X[row_x][opp_move_col] = opponent

            can_O_win_now = False
            for move2_col in get_valid_moves(board_after_X):
                board_after_O2 = copy.deepcopy(board_after_X)
                row2 = get_next_open_row(board_after_O2, move2_col)
                board_after_O2[row2][move2_col] = player
                if check_win(board_after_O2, player):
                    can_O_win_now = True
                    break
            
            if not can_O_win_now:
                is_forced_win = False
                break
        
        if is_forced_win:
            win_in_2_moves.append(f"{col_to_char(move_col)}{row_to_char(row1)}")

    if win_in_2_moves:
        print(', '.join(sorted(win_in_2_moves)))
    else:
        # Fallback if no quick wins are found
        print("No optimal winning moves found in 1 or 2 turns.")


solve_connect4()