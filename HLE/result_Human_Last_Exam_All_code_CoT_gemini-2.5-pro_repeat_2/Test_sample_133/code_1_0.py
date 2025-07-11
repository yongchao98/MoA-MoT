def solve_connect_four():
    """
    This function analyzes the given Connect 4 board state to find the optimal
    moves for player 'O' to win as fast as possible.
    """
    # The board is represented as a list of lists.
    # Row 0 is the top row as shown in the problem's image.
    board = [
        ['.', '.', '.', '.', '.', '.', '.'],  # Row 1
        ['.', '.', '.', '.', '.', '.', '.'],  # Row 2
        ['.', '.', '.', '.', '.', '.', '.'],  # Row 3
        ['.', '.', '.', 'O', 'O', '.', '.'],  # Row 4
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],  # Row 5
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']   # Row 6
    ]
    
    player = 'O'
    opponent = 'X'
    rows, cols = 6, 7

    def check_win(b, p):
        # Check horizontal
        for c in range(cols - 3):
            for r in range(rows):
                if b[r][c] == p and b[r][c+1] == p and b[r][c+2] == p and b[r][c+3] == p:
                    return True
        # Check vertical
        for c in range(cols):
            for r in range(rows - 3):
                if b[r][c] == p and b[r+1][c] == p and b[r+2][c] == p and b[r+3][c] == p:
                    return True
        # Check positively sloped diagonals
        for c in range(cols - 3):
            for r in range(rows - 3):
                if b[r][c] == p and b[r+1][c+1] == p and b[r+2][c+2] == p and b[r+3][c+3] == p:
                    return True
        # Check negatively sloped diagonals
        for c in range(cols - 3):
            for r in range(3, rows):
                if b[r][c] == p and b[r-1][c+1] == p and b[r-2][c+2] == p and b[r-3][c+3] == p:
                    return True
        return False

    def get_next_open_row(b, c):
        for r in range(rows - 1, -1, -1):
            if b[r][c] == '.':
                return r
        return None

    def get_valid_cols(b):
        return [c for c in range(cols) if b[0][c] == '.']

    # Tier 1: Check for win in 1 move
    winning_moves_1_ply = []
    for col in get_valid_cols(board):
        row = get_next_open_row(board, col)
        temp_board = [r[:] for r in board]
        temp_board[row][col] = player
        if check_win(temp_board, player):
            move_str = f"{chr(ord('a') + col)}{row + 1}"
            winning_moves_1_ply.append(move_str)
    
    if winning_moves_1_ply:
        print(", ".join(sorted(winning_moves_1_ply)))
        return

    # Tier 2: Check for win in 3 moves (forced win)
    winning_moves_3_ply = []
    for col in get_valid_cols(board):
        row = get_next_open_row(board, col)
        
        board_after_o1 = [r[:] for r in board]
        board_after_o1[row][col] = player

        # Skip this move if it allows the opponent to win immediately
        can_opponent_win = False
        for opp_col in get_valid_cols(board_after_o1):
             opp_row = get_next_open_row(board_after_o1, opp_col)
             board_after_x = [r[:] for r in board_after_o1]
             board_after_x[opp_row][opp_col] = opponent
             if check_win(board_after_x, opponent):
                 can_opponent_win = True
                 break
        if can_opponent_win:
            continue
            
        # Count winning threats created by this move. A threat is an empty
        # spot where O can win on the next turn. We need >= 2 threats.
        threat_cols = set()
        for win_col in get_valid_cols(board_after_o1):
            win_row = get_next_open_row(board_after_o1, win_col)
            
            win_board = [r[:] for r in board_after_o1]
            win_board[win_row][win_col] = player
            if check_win(win_board, player):
                threat_cols.add(win_col)
        
        if len(threat_cols) >= 2:
            # Found a move that creates an unblockable win.
            # Convert to "a1" notation. Row numbering is 1-indexed from the top.
            move_str = f"{chr(ord('a') + col)}{row + 1}"
            winning_moves_3_ply.append(move_str)
            
    if winning_moves_3_ply:
        print(", ".join(sorted(winning_moves_3_ply)))
        return
        
    print("No optimal winning moves found.")

solve_connect_four()