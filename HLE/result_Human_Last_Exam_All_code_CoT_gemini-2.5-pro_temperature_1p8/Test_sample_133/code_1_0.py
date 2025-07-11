import numpy as np

# Define constants
ROW_COUNT = 6
COLUMN_COUNT = 7
PLAYER_O = 'O'
PLAYER_X = 'X'
EMPTY = '.'
WINDOW_LENGTH = 4

def create_board():
    """Creates the initial board state as specified in the problem."""
    board = np.full((ROW_COUNT, COLUMN_COUNT), EMPTY, dtype=str)
    # Using the coordinate system where (0,0) is top-left
    # From image: 4d=O, 4e=O
    board[3][3] = PLAYER_O
    board[3][4] = PLAYER_O
    # 5a=O, 5c=X, 5d=O, 5e=X, 5f=X, 5g=X
    board[4][0] = PLAYER_O
    board[4][2] = PLAYER_X
    board[4][3] = PLAYER_O
    board[4][4] = PLAYER_X
    board[4][5] = PLAYER_X
    board[4][6] = PLAYER_X
    # 6a=X, 6b=O, 6c=O, 6d=X, 6e=X, 6f=O, 6g=X
    board[5][0] = PLAYER_X
    board[5][1] = PLAYER_O
    board[5][2] = PLAYER_O
    board[5][3] = PLAYER_X
    board[5][4] = PLAYER_X
    board[5][5] = PLAYER_O
    board[5][6] = PLAYER_X
    return board

def is_valid_location(board, col):
    """Checks if a column is not full."""
    return board[0][col] == EMPTY

def get_next_open_row(board, col):
    """Returns the row index for the next piece in a column."""
    for r in range(ROW_COUNT - 1, -1, -1):
        if board[r][col] == EMPTY:
            return r
    return None

def drop_piece(board, row, col, piece):
    """Places a piece on the board."""
    board[row][col] = piece

def check_win(board, piece):
    """Checks if a player has won."""
    # Check horizontal
    for r in range(ROW_COUNT):
        for c in range(COLUMN_COUNT - 3):
            if all(board[r][c+i] == piece for i in range(WINDOW_LENGTH)):
                return True
    # Check vertical
    for c in range(COLUMN_COUNT):
        for r in range(ROW_COUNT - 3):
            if all(board[r+i][c] == piece for i in range(WINDOW_LENGTH)):
                return True
    # Check positive diagonal
    for c in range(COLUMN_COUNT - 3):
        for r in range(ROW_COUNT - 3):
            if all(board[r+i][c+i] == piece for i in range(WINDOW_LENGTH)):
                return True
    # Check negative diagonal
    for c in range(COLUMN_COUNT - 3):
        for r in range(3, ROW_COUNT):
            if all(board[r-i][c+i] == piece for i in range(WINDOW_LENGTH)):
                return True
    return False

def get_move_name(col, row):
    """Converts (col, row) indices to human-readable format like 'a1'."""
    col_name = chr(ord('a') + col)
    # The visual grid row is 1-indexed, our array is 0-indexed.
    row_name = str(row + 1)
    return f"{col_name}{row_name}"

def solve():
    """Finds the optimal moves for player 'O'."""
    board = create_board()
    player = PLAYER_O

    # --- Step 1: Check for any immediate winning moves ---
    immediate_wins = []
    for col in range(COLUMN_COUNT):
        if is_valid_location(board, col):
            row = get_next_open_row(board, col)
            board_copy = board.copy()
            drop_piece(board_copy, row, col, player)
            if check_win(board_copy, player):
                immediate_wins.append(get_move_name(col, row))
    
    if immediate_wins:
        print(", ".join(sorted(immediate_wins)))
        return

    # --- Step 2: Look for a forced win (creating multiple unblockable threats) ---
    forced_wins = []
    for col1 in range(COLUMN_COUNT):
        if is_valid_location(board, col1):
            row1 = get_next_open_row(board, col1)
            board_after_o1 = board.copy()
            drop_piece(board_after_o1, row1, col1, player)

            # Find how many winning threats are created by this move.
            # A threat is a spot where O can play next turn to win.
            winning_threat_cols = set()
            for col2 in range(COLUMN_COUNT):
                if is_valid_location(board_after_o1, col2):
                    row2 = get_next_open_row(board_after_o1, col2)
                    board_after_o2 = board_after_o1.copy()
                    
                    drop_piece(board_after_o2, row2, col2, player)
                    if check_win(board_after_o2, player):
                        winning_threat_cols.add(col2)
            
            # If threats exist in at least two different columns, the opponent
            # cannot block both. This move is a forced win.
            if len(winning_threat_cols) >= 2:
                forced_wins.append(get_move_name(col1, row1))

    if forced_wins:
        print(", ".join(sorted(forced_wins)))
        return
        
    print("No optimal winning move found with the current strategy.")

if __name__ == '__main__':
    solve()