import sys

# Constants
ROWS = 6
COLS = 7
PLAYER_O = 'O'
PLAYER_X = 'X'
EMPTY = '.'

def create_board():
    """Initializes the game board based on the problem description."""
    board = [[EMPTY for _ in range(COLS)] for _ in range(ROWS)]
    board[3][3] = PLAYER_O # d4
    board[3][4] = PLAYER_O # e4
    board[4][0] = PLAYER_O # a5
    board[4][2] = PLAYER_X # c5
    board[4][3] = PLAYER_O # d5
    board[4][4] = PLAYER_X # e5
    board[4][5] = PLAYER_X # f5
    board[4][6] = PLAYER_X # g5
    board[5][0] = PLAYER_X # a6
    board[5][1] = PLAYER_O # b6
    board[5][2] = PLAYER_O # c6
    board[5][3] = PLAYER_X # d6
    board[5][4] = PLAYER_X # e6
    board[5][5] = PLAYER_O # f6
    board[5][6] = PLAYER_X # g6
    return board

def get_next_open_row(board, col):
    """Returns the lowest empty row number for a given column."""
    for r in range(ROWS - 1, -1, -1):
        if board[r][col] == EMPTY:
            return r
    return None

def is_winning_move(board, player, r, c):
    """Checks if placing a piece at (r, c) results in a win."""
    # Check horizontal
    for c_start in range(COLS - 3):
        if all(board[r][c_start + i] == player for i in range(4)):
            return True
            
    # Check vertical
    for r_start in range(ROWS - 3):
        if all(board[r_start + i][c] == player for i in range(4)):
            return True

    # Check positive diagonal (\)
    for r_start in range(ROWS - 3):
        for c_start in range(COLS - 3):
            if all(board[r_start + i][c_start + i] == player for i in range(4)):
                return True

    # Check negative diagonal (/)
    for r_start in range(3, ROWS):
        for c_start in range(COLS - 3):
            if all(board[r_start - i][c_start + i] == player for i in range(4)):
                return True
                
    return False

def find_optimal_moves():
    """Finds all optimal moves for the current player."""
    board = create_board()
    optimal_moves = []
    
    # First, check for any immediate winning moves (win in 1)
    for col in range(COLS):
        row = get_next_open_row(board, col)
        if row is not None:
            board[row][col] = PLAYER_O
            if is_winning_move(board, PLAYER_O, row, col):
                # An immediate win is the fastest possible
                move_name = f"{chr(ord('a') + col)}{ROWS - row}"
                optimal_moves.append(move_name)
            board[row][col] = EMPTY # backtrack
    
    if optimal_moves:
        print(", ".join(sorted(optimal_moves)))
        return

    # If no immediate win, check for forced wins (win in 3)
    # A move is optimal if it creates at least two unstoppable threats for the next turn.
    for col in range(COLS):
        row = get_next_open_row(board, col)
        if row is not None:
            # Simulate O's move
            board[row][col] = PLAYER_O
            
            # Count how many winning moves O has on the next turn
            winning_threats = 0
            for next_col in range(COLS):
                next_row = get_next_open_row(board, next_col)
                if next_row is not None:
                    # Check if placing a piece here would win
                    board[next_row][next_col] = PLAYER_O
                    if is_winning_move(board, PLAYER_O, next_row, next_col):
                        winning_threats += 1
                    board[next_row][next_col] = EMPTY # backtrack check
            
            # If there are 2 or more threats, X cannot block them all.
            if winning_threats >= 2:
                move_name = f"{chr(ord('a') + col)}{ROWS - row}"
                optimal_moves.append(move_name)

            # Backtrack O's initial move
            board[row][col] = EMPTY
            
    print(", ".join(sorted(optimal_moves)))

find_optimal_moves()