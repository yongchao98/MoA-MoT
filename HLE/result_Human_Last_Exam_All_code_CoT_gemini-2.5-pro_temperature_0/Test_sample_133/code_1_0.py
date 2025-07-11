import numpy as np

# Constants
ROWS = 6
COLS = 7
PLAYER_O = 'O'
PLAYER_X = 'X'
EMPTY = '.'

def create_board():
    """Creates the initial board state for the problem."""
    board = np.full((ROWS, COLS), EMPTY)
    # The board is represented with row 0 at the top and row 5 at the bottom.
    # The user's board is 1-indexed from the top.
    # Row 6 in the prompt corresponds to index 5 here.
    board[5] = ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    # Row 5 in the prompt corresponds to index 4 here.
    board[4] = ['O', '.', 'X', 'O', 'X', 'X', 'X']
    # Row 4 in the prompt corresponds to index 3 here.
    board[3] = ['.', '.', '.', 'O', 'O', '.', '.']
    return board

def is_valid_location(board, col):
    """Checks if a column is not full."""
    return board[0][col] == EMPTY

def get_next_open_row(board, col):
    """Gets the row a piece will fall to in a column."""
    for r in range(ROWS - 1, -1, -1):
        if board[r][col] == EMPTY:
            return r
    return None

def drop_piece(board, row, col, piece):
    """Drops a piece on the board."""
    board[row][col] = piece

def winning_move(board, piece):
    """Checks if the given piece has a winning configuration."""
    # Check horizontal locations for win
    for c in range(COLS - 3):
        for r in range(ROWS):
            if all(board[r][c+i] == piece for i in range(4)):
                return True

    # Check vertical locations for win
    for c in range(COLS):
        for r in range(ROWS - 3):
            if all(board[r+i][c] == piece for i in range(4)):
                return True

    # Check positively sloped diagonals
    for c in range(COLS - 3):
        for r in range(ROWS - 3):
            if all(board[r+i][c+i] == piece for i in range(4)):
                return True

    # Check negatively sloped diagonals
    for c in range(COLS - 3):
        for r in range(3, ROWS):
            if all(board[r-i][c+i] == piece for i in range(4)):
                return True
    return False

def get_char_from_col(col):
    """Converts a column index (0-6) to a character ('a'-'g')."""
    return chr(ord('a') + col)

def solve_connect4():
    """
    Finds and prints the optimal moves for 'O' to win as fast as possible.
    """
    initial_board = create_board()
    
    # --- Level 1: Check for win in 1 move ---
    win_in_1_moves = []
    for col in range(COLS):
        if is_valid_location(initial_board, col):
            board_copy = initial_board.copy()
            row = get_next_open_row(board_copy, col)
            drop_piece(board_copy, row, col, PLAYER_O)
            if winning_move(board_copy, PLAYER_O):
                # Convert 0-indexed row to 1-indexed row from top
                win_in_1_moves.append(f"{get_char_from_col(col)}{row + 1}")

    if win_in_1_moves:
        print(", ".join(win_in_1_moves))
        return

    # --- Level 2: Check for win in 2 moves (forced win) ---
    win_in_2_moves = []
    for col_o1 in range(COLS):
        if not is_valid_location(initial_board, col_o1):
            continue

        board_after_o1 = initial_board.copy()
        row_o1 = get_next_open_row(board_after_o1, col_o1)
        drop_piece(board_after_o1, row_o1, col_o1, PLAYER_O)

        # A move is a forced win if for ALL of X's possible responses, O has a winning move.
        is_forced_win = True
        
        # Check if there are any valid moves for X.
        x_has_valid_moves = any(is_valid_location(board_after_o1, c) for c in range(COLS))
        if not x_has_valid_moves:
            is_forced_win = False 

        for col_x1 in range(COLS):
            if not is_valid_location(board_after_o1, col_x1):
                continue

            board_after_x1 = board_after_o1.copy()
            row_x1 = get_next_open_row(board_after_x1, col_x1)
            drop_piece(board_after_x1, row_x1, col_x1, PLAYER_X)
            
            # If X can win, then O's initial move was a blunder.
            if winning_move(board_after_x1, PLAYER_X):
                is_forced_win = False
                break

            # After X moves, does O have a guaranteed winning response?
            o_can_win_in_response = False
            for col_o2 in range(COLS):
                if is_valid_location(board_after_x1, col_o2):
                    board_after_o2 = board_after_x1.copy()
                    row_o2 = get_next_open_row(board_after_o2, col_o2)
                    drop_piece(board_after_o2, row_o2, col_o2, PLAYER_O)
                    if winning_move(board_after_o2, PLAYER_O):
                        o_can_win_in_response = True
                        break
            
            if not o_can_win_in_response:
                is_forced_win = False
                break
        
        if is_forced_win:
            # Convert 0-indexed row to 1-indexed row from top for user output
            win_in_2_moves.append(f"{get_char_from_col(col_o1)}{ROWS - row_o1}")

    if win_in_2_moves:
        print(", ".join(win_in_2_moves))
        return
        
    print("No optimal winning move found within 2 turns.")

if __name__ == '__main__':
    solve_connect4()