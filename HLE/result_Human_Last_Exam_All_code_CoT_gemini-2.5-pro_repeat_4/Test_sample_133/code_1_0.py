import copy

def solve_connect_four():
    """
    Analyzes the given Connect 4 board to find all optimal moves for 'O'.
    An optimal move is one that forces a win in the minimum number of turns.
    """

    # Board representation: 6 rows, 7 columns. (0,0) is bottom-left.
    # Diagram rows 1-6 (top-down) map to my rows 5-0 (top-down).
    # Diagram columns a-g map to my columns 0-6.
    board = [
        ['X', 'O', 'O', 'X', 'X', 'O', 'X'],  # Row 6 in diagram
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],  # Row 5 in diagram
        ['.', '.', '.', 'O', 'O', '.', '.'],  # Row 4 in diagram
        ['.', '.', '.', '.', '.', '.', '.'],  # Row 3 in diagram
        ['.', '.', '.', '.', '.', '.', '.'],  # Row 2 in diagram
        ['.', '.', '.', '.', '.', '.', '.']   # Row 1 in diagram
    ]
    ROWS = 6
    COLS = 7
    PLAYER = 'O'
    OPPONENT = 'X'

    def get_next_open_row(board, col):
        """Finds the first empty row in a given column."""
        for r in range(ROWS):
            if board[r][col] == '.':
                return r
        return -1 # Column is full

    def get_valid_moves(board):
        """Returns a list of all valid moves as (row, col) tuples."""
        valid_moves = []
        for c in range(COLS):
            r = get_next_open_row(board, c)
            if r != -1:
                valid_moves.append((r, c))
        return valid_moves

    def check_win(board, player, r, c):
        """Checks if the move at (r, c) by player results in a win."""
        # Check horizontal
        for i in range(COLS - 3):
            if board[r][i] == player and board[r][i+1] == player and board[r][i+2] == player and board[r][i+3] == player:
                return True
        # Check vertical
        for i in range(ROWS - 3):
            if board[i][c] == player and board[i+1][c] == player and board[i+2][c] == player and board[i+3][c] == player:
                return True
        # Check positive diagonal
        for i in range(ROWS - 3):
            for j in range(COLS - 3):
                if board[i][j] == player and board[i+1][j+1] == player and board[i+2][j+2] == player and board[i+3][j+3] == player:
                    return True
        # Check negative diagonal
        for i in range(3, ROWS):
            for j in range(COLS - 3):
                if board[i][j] == player and board[i-1][j+1] == player and board[i-2][j+2] == player and board[i-3][j+3] == player:
                    return True
        return False

    def find_threats(board, player):
        """Finds all winning moves (threats) for the player."""
        threats = set()
        for r, c in get_valid_moves(board):
            temp_board = copy.deepcopy(board)
            temp_board[r][c] = player
            if check_win(temp_board, player, r, c):
                threats.add((r, c))
        return threats

    optimal_moves_coords = []
    
    # There are no immediate winning moves, so we search for forced wins in 2 turns.
    initial_valid_moves = get_valid_moves(board)

    for r1, c1 in initial_valid_moves:
        # 1. O makes a move
        board_after_O1 = copy.deepcopy(board)
        board_after_O1[r1][c1] = PLAYER
        
        # Find O's threats after this move
        threats_for_O = find_threats(board_after_O1, PLAYER)
        
        if not threats_for_O:
            continue

        # 2. Assume X will try to block. A move is optimal for O if for ALL of X's
        # possible counter-moves, O still has a winning move left.
        is_defended = False
        opponent_valid_moves = get_valid_moves(board_after_O1)
        
        if not opponent_valid_moves: # No more moves left for opponent
             is_defended = True

        for r2, c2 in opponent_valid_moves:
            # X makes a counter-move
            board_after_X1 = copy.deepcopy(board_after_O1)
            board_after_X1[r2][c2] = OPPONENT
            
            # If X's move leaves O with NO winning threats, then X has defended.
            if not find_threats(board_after_X1, PLAYER):
                is_defended = True
                break
        
        if not is_defended:
            optimal_moves_coords.append((r1, c1))

    # Convert coordinates to algebraic notation (e.g., (2,2) -> c4)
    def to_alg_notation(r, c):
        col_char = chr(ord('a') + c)
        row_char = str(6 - r)
        return f"{col_char}{row_char}"

    result = sorted([to_alg_notation(r,c) for r,c in optimal_moves_coords])
    print(", ".join(result))

solve_connect_four()