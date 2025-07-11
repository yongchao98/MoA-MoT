def solve_connect4():
    """
    Analyzes the given Connect 4 board to find the optimal move for 'O'.
    The function checks for immediate wins and then for forced wins (next turn).
    """

    # Board representation: 0,0 is bottom-left
    # '.' is empty, 'X' and 'O' are players.
    # The board is represented with rows 0-5 (bottom to top)
    # and columns 0-6 (left to right).
    board = [
      ['X','O','O','X','X','O','X'],  # row 0 (board's row 6)
      ['O','.','X','O','X','X','X'],  # row 1 (board's row 5)
      ['.','.','.','O','O','.','.'],  # row 2 (board's row 4)
      ['.','.','.','.','.','.','.'],  # row 3 (board's row 3)
      ['.','.','.','.','.','.','.'],  # row 4 (board's row 2)
      ['.','.','.','.','.','.','.']   # row 5 (board's row 1)
    ]
    ROWS = 6
    COLS = 7
    PLAYER = 'O'
    OPPONENT = 'X'

    def get_valid_moves(current_board):
        """Finds the next available row for each non-full column."""
        moves = []
        for c in range(COLS):
            for r in range(ROWS):
                if current_board[r][c] == '.':
                    moves.append((r, c))
                    break
        return moves

    def check_win(current_board, player, r, c):
        """Checks if placing a piece at (r, c) results in a win."""
        # Horizontal
        count = 0
        for i in range(COLS):
            if current_board[r][i] == player:
                count += 1
            else:
                count = 0
            if count >= 4: return True
        # Vertical
        count = 0
        for i in range(ROWS):
            if current_board[i][c] == player:
                count += 1
            else:
                count = 0
            if count >= 4: return True
        # Diagonal /
        count = 0
        for i in range(-3, 4):
            if 0 <= r + i < ROWS and 0 <= c + i < COLS:
                if current_board[r + i][c + i] == player:
                    count += 1
                else:
                    count = 0
                if count >= 4: return True
        # Diagonal \
        count = 0
        for i in range(-3, 4):
            if 0 <= r + i < ROWS and 0 <= c - i < COLS:
                if current_board[r + i][c - i] == player:
                    count += 1
                else:
                    count = 0
                if count >= 4: return True
        return False

    # Step 1: Check for immediate winning moves
    winning_moves = []
    valid_moves = get_valid_moves(board)
    for r, c in valid_moves:
        temp_board = [row[:] for row in board]
        temp_board[r][c] = PLAYER
        if check_win(temp_board, PLAYER, r, c):
            winning_moves.append((r, c))

    if winning_moves:
        # Convert to problem's notation (e.g., f4)
        results = []
        for r, c in winning_moves:
            col_char = chr(ord('a') + c)
            row_num = 6 - r
            results.append(f"{col_char}{row_num}")
        print(", ".join(results))
        return

    # Step 2: Check for forced wins on the next turn
    forced_win_moves = []
    for r1, c1 in valid_moves: # O's first move
        board1 = [row[:] for row in board]
        board1[r1][c1] = PLAYER
        
        # Assume this move leads to a win, unless opponent can prevent it
        can_be_countered = False
        
        opponent_moves = get_valid_moves(board1)
        if not opponent_moves: # Draw
             continue

        for r2, c2 in opponent_moves: # X's response
            board2 = [row[:] for row in board1]
            board2[r2][c2] = OPPONENT
            
            # Can O win on the next turn, regardless of X's move?
            o_can_win_now = False
            o_winning_responses = get_valid_moves(board2)
            for r3, c3 in o_winning_responses:
                board3 = [row[:] for row in board2]
                board3[r3][c3] = PLAYER
                if check_win(board3, PLAYER, r3, c3):
                    o_can_win_now = True
                    break
            
            if not o_can_win_now:
                # X found a move to stop an immediate win from O
                can_be_countered = True
                break
        
        if not can_be_countered:
            forced_win_moves.append((r1, c1))

    if forced_win_moves:
        results = []
        for r, c in forced_win_moves:
            col_char = chr(ord('a') + c)
            row_num = 6 - r
            results.append(f"{col_char}{row_num}")
        print(", ".join(results))
        return

    print("No optimal move to win found.")

solve_connect4()