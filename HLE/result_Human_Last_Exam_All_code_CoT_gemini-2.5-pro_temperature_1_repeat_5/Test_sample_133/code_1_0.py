import copy

# Constants for the board representation
ROWS = 6
COLS = 7
EMPTY = 0
X_PIECE = 1
O_PIECE = 2

def create_board():
    """Creates the initial board state from the problem description."""
    # Mapping: R1->row 0, R6->row 5; a->col 0, g->col 6
    board = [
        [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY], # R1
        [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY], # R2
        [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY], # R3
        [EMPTY, EMPTY, EMPTY, O_PIECE, O_PIECE, EMPTY, EMPTY], # R4
        [O_PIECE, EMPTY, X_PIECE, O_PIECE, X_PIECE, X_PIECE, X_PIECE], # R5
        [X_PIECE, O_PIECE, O_PIECE, X_PIECE, X_PIECE, O_PIECE, X_PIECE]  # R6
    ]
    return board

def get_next_open_row(board, col):
    """Returns the next open row in a given column."""
    for r in range(ROWS - 1, -1, -1):
        if board[r][col] == EMPTY:
            return r
    return None

def drop_piece(board, row, col, piece):
    """Drops a piece in the board and returns a new board."""
    board_copy = copy.deepcopy(board)
    board_copy[row][col] = piece
    return board_copy

def check_winning_move(board, piece):
    """Checks if the specified player has won."""
    # Check horizontal
    for c in range(COLS - 3):
        for r in range(ROWS):
            if board[r][c] == piece and board[r][c+1] == piece and board[r][c+2] == piece and board[r][c+3] == piece:
                return True
    # Check vertical
    for c in range(COLS):
        for r in range(ROWS - 3):
            if board[r][c] == piece and board[r+1][c] == piece and board[r+2][c] == piece and board[r+3][c] == piece:
                return True
    # Check positively sloped diagonals
    for c in range(COLS - 3):
        for r in range(ROWS - 3):
            if board[r][c] == piece and board[r+1][c+1] == piece and board[r+2][c+2] == piece and board[r+3][c+3] == piece:
                return True
    # Check negatively sloped diagonals
    for c in range(COLS - 3):
        for r in range(3, ROWS):
            if board[r][c] == piece and board[r-1][c+1] == piece and board[r-2][c+2] == piece and board[r-3][c+3] == piece:
                return True
    return False

def get_valid_cols(board):
    """Gets a list of columns that are not full."""
    valid_cols = []
    for col in range(COLS):
        if board[0][col] == EMPTY:
            valid_cols.append(col)
    return valid_cols

def col_to_letter(col):
    """Converts a column index to its letter representation."""
    return chr(ord('a') + col)

def find_optimal_moves():
    """Analyzes the board to find all optimal moves for 'O' to win as fast as possible."""
    initial_board = create_board()
    player = O_PIECE
    opponent = X_PIECE
    optimal_moves = []

    possible_moves = get_valid_cols(initial_board)

    for col in possible_moves:
        row = get_next_open_row(initial_board, col)
        if row is None:
            continue

        # Simulate O's move
        board_after_O1 = drop_piece(initial_board, row, col, player)

        # Win in 2 moves: Check if for every possible opponent move, we have a winning reply.
        can_opponent_defend = False
        opponent_moves = get_valid_cols(board_after_O1)
        
        # If there are no more moves for the opponent, it's a draw, not a win for O
        if not opponent_moves:
            continue

        for opp_col in opponent_moves:
            opp_row = get_next_open_row(board_after_O1, opp_col)
            board_after_X = drop_piece(board_after_O1, opp_row, opp_col, opponent)

            # Does O have a winning reply?
            has_winning_reply = False
            our_next_moves = get_valid_cols(board_after_X)
            for our_next_col in our_next_moves:
                our_next_row = get_next_open_row(board_after_X, our_next_col)
                final_board = drop_piece(board_after_X, our_next_row, our_next_col, player)
                if check_winning_move(final_board, player):
                    has_winning_reply = True
                    break # Found a winning reply for this opponent move

            # If for this one opponent move, we have no winning reply, this path is not a guaranteed win.
            if not has_winning_reply:
                can_opponent_defend = True
                break # Opponent can defend, so this initial move is not optimal (win in 2)

        if not can_opponent_defend:
            # This move leads to a win regardless of what the opponent does.
            move_name = f"{col_to_letter(col)}{ROWS - row}"
            optimal_moves.append(move_name)

    print(", ".join(sorted(optimal_moves)))

if __name__ == "__main__":
    find_optimal_moves()