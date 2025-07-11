import copy

def find_optimal_connect4_moves():
    """
    Analyzes a specific Connect 4 board to find the optimal moves for player 'O'.
    An optimal move is one that secures a win in the minimum number of turns.
    """
    # Board representation: rows 0-5 (top to bottom), cols 0-6 (left to right)
    # The image's Row 1 is index 0, Row 6 is index 5.
    # Column 'a' is index 0, 'g' is index 6.
    initial_board = [
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', 'O', '.', '.'],
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    ]
    player_o = 'O'
    player_x = 'X'
    rows, cols = 6, 7

    def get_next_available_row(board, col):
        for r in range(rows - 1, -1, -1):
            if board[r][col] == '.':
                return r
        return -1

    def make_move(board, col, player):
        temp_board = copy.deepcopy(board)
        row = get_next_available_row(temp_board, col)
        if row != -1:
            temp_board[row][col] = player
        return temp_board

    def check_win(board, player):
        # Check horizontal
        for r in range(rows):
            for c in range(cols - 3):
                if all(board[r][c+i] == player for i in range(4)):
                    return True
        # Check vertical
        for c in range(cols):
            for r in range(rows - 3):
                if all(board[r+i][c] == player for i in range(4)):
                    return True
        # Check diagonal / (bottom-left to top-right)
        for r in range(3, rows):
            for c in range(cols - 3):
                if all(board[r-i][c+i] == player for i in range(4)):
                    return True
        # Check diagonal \ (top-left to bottom-right)
        for r in range(rows - 3):
            for c in range(cols - 3):
                if all(board[r+i][c+i] == player for i in range(4)):
                    return True
        return False

    def get_legal_moves(board):
        return [c for c in range(cols) if get_next_available_row(board, c) != -1]

    def col_to_letter(col):
        return chr(ord('a') + col)

    # Step 1: Check for immediate wins (win in 1 move).
    # Based on analysis, none exist, but we check to be thorough.
    winning_moves_in_1 = []
    for move_col in get_legal_moves(initial_board):
        board_after_move = make_move(initial_board, move_col, player_o)
        if check_win(board_after_move, player_o):
            winning_moves_in_1.append(col_to_letter(move_col))
    
    if winning_moves_in_1:
        print(", ".join(sorted(winning_moves_in_1)))
        return

    # Step 2: Check for forced wins in 3 moves (O -> X -> O win).
    winning_moves_in_3 = []
    for move_o1_col in get_legal_moves(initial_board):
        board_after_o1 = make_move(initial_board, move_o1_col, player_o)
        
        is_forcing_move = True
        opponent_legal_moves = get_legal_moves(board_after_o1)

        if not opponent_legal_moves:
            continue

        for move_x_col in opponent_legal_moves:
            board_after_x = make_move(board_after_o1, move_x_col, player_x)
            
            o_can_win = False
            for move_o2_col in get_legal_moves(board_after_x):
                board_after_o2 = make_move(board_after_x, move_o2_col, player_o)
                if check_win(board_after_o2, player_o):
                    o_can_win = True
                    break
            
            if not o_can_win:
                is_forcing_move = False
                break
        
        if is_forcing_move:
            winning_moves_in_3.append(col_to_letter(move_o1_col))

    print(", ".join(sorted(winning_moves_in_3)))

# Execute the function to find and print the optimal moves.
find_optimal_connect4_moves()