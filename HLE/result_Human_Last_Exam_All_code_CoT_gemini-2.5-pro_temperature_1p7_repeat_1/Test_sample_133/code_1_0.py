import copy

def solve_connect4():
    """
    Analyzes a Connect 4 board to find the optimal moves for player 'O'.
    """
    # Board setup: 6 rows, 7 columns. '.' is empty, 'X' and 'O' are players.
    # (0,0) is top-left (a1).
    board = [
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', 'O', '.', '.'],
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    ]
    PLAYER = 'O'
    OPPONENT = 'X'
    NUM_ROWS = 6
    NUM_COLS = 7

    def get_landing_row(board_state, col):
        """Find the lowest empty row in a given column."""
        for r in range(NUM_ROWS - 1, -1, -1):
            if board_state[r][col] == '.':
                return r
        return -1 # Column is full

    def check_win_at(board_state, player, r, c):
        """Check if placing a piece at (r, c) results in a win."""
        # Check horizontal
        for col_start in range(max(0, c - 3), min(NUM_COLS - 3, c + 1)):
            if all(board_state[r][col_start + i] == player for i in range(4)):
                return True
        # Check vertical
        for row_start in range(max(0, r - 3), min(NUM_ROWS - 3, r + 1)):
            if all(board_state[row_start + i][c] == player for i in range(4)):
                return True
        # Check diagonal /
        for i in range(4):
            start_r, start_c = r + i, c - i
            if 0 <= start_r <= NUM_ROWS - 4 and 0 <= start_c <= NUM_COLS - 4:
                if all(board_state[start_r - j][start_c + j] == player for j in range(4)):
                    return True
        # Check diagonal \
        for i in range(4):
            start_r, start_c = r - i, c - i
            if 0 <= start_r <= NUM_ROWS - 4 and 0 <= start_c <= NUM_COLS - 4:
                 if all(board_state[start_r + j][start_c + j] == player for j in range(4)):
                    return True
        return False

    def to_notation(r, c):
        """Converts (row, col) coordinates to game notation like 'a1'."""
        col_char = chr(ord('a') + c)
        row_char = str(r + 1)
        return f"{col_char}{row_char}"

    # --- Main Logic ---

    # Level 1: Find any moves that win in 1 turn.
    winning_moves_now = []
    for c in range(NUM_COLS):
        r = get_landing_row(board, c)
        if r != -1:
            temp_board = copy.deepcopy(board)
            temp_board[r][c] = PLAYER
            if check_win_at(temp_board, PLAYER, r, c):
                winning_moves_now.append(to_notation(r,c))
    
    if winning_moves_now:
        print(", ".join(sorted(winning_moves_now)))
        return

    # Level 2: Find any moves that guarantee a win in 3 turns (forks).
    guaranteed_wins = []
    for c in range(NUM_COLS):
        r = get_landing_row(board, c)
        if r == -1:
            continue

        # Simulate O's move
        board_after_O = copy.deepcopy(board)
        board_after_O[r][c] = PLAYER

        # Can opponent win immediately after our move? If so, it's a bad move.
        can_opponent_win = False
        for c_opp in range(NUM_COLS):
            r_opp = get_landing_row(board_after_O, c_opp)
            if r_opp != -1:
                board_after_opp_move = copy.deepcopy(board_after_O)
                board_after_opp_move[r_opp][c_opp] = OPPONENT
                if check_win_at(board_after_opp_move, OPPONENT, r_opp, c_opp):
                    can_opponent_win = True
                    break
        if can_opponent_win:
             continue # This move allows opponent to win, so it's not optimal.


        # Find all of O's winning spots on the next turn
        winning_spots_for_O = set()
        for c_win in range(NUM_COLS):
            r_win = get_landing_row(board_after_O, c_win)
            if r_win != -1:
                board_to_check_win = copy.deepcopy(board_after_O)
                board_to_check_win[r_win][c_win] = PLAYER
                if check_win_at(board_to_check_win, PLAYER, r_win, c_win):
                    winning_spots_for_O.add((r_win, c_win))

        # If there are at least two threats that opponent cannot block in one move
        if len(winning_spots_for_O) >= 2:
            # Check if threats are in different columns, which is a standard fork
            # A more complex check would see if one move by X could block multiple threats
            # (e.g., if threats are vertically aligned). Here, we check for a simpler fork.
            unblockable_threat_count = 0
            
            # Simulate X's blocking attempts
            is_guaranteed = True
            for c_block in range(NUM_COLS):
                r_block = get_landing_row(board_after_O, c_block)
                if r_block == -1:
                    continue
                
                # Check remaining threats after X plays at (r_block, c_block)
                remaining_threats = winning_spots_for_O.copy()
                if (r_block, c_block) in remaining_threats:
                    remaining_threats.remove((r_block, c_block))

                # What if X's move opens up a new winning spot for O?
                # This is complex, but the current identified threats should suffice.
                # Here, we assume X's move only serves to block.
                
                if not remaining_threats:
                    is_guaranteed = False
                    break
            
            if is_guaranteed and winning_spots_for_O:
                guaranteed_wins.append(to_notation(r, c))

    print(", ".join(sorted(guaranteed_wins)))

solve_connect4()
>>> f4