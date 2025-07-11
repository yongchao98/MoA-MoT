import copy

def solve_connect_four():
    """
    Analyzes a Connect 4 board to find the optimal winning moves for player 'O'.
    """
    # Board setup: row 0 is the top, row 5 is the bottom
    # The input board is visually flipped, so we map it carefully.
    # [ . . . . . . . ] (row 0)
    # [ . . . . . . . ] (row 1)
    # [ . . . . . . . ] (row 2)
    # [ . . . O O . . ] (row 3) -> corresponds to row 4 in the puzzle image
    # [ O . X O X X X ] (row 4) -> corresponds to row 5 in the puzzle image
    # [ X O O X X O X ] (row 5) -> corresponds to row 6 in the puzzle image
    #   a b c d e f g
    
    board = [
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', 'O', '.', '.'],
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    ]

    ROWS = 6
    COLS = 7
    ME = 'O'
    OPPONENT = 'X'

    def get_next_open_row(b, col):
        """Finds the first empty row in a given column."""
        for r in range(ROWS - 1, -1, -1):
            if b[r][col] == '.':
                return r
        return None

    def check_win(b, player):
        """Checks if the given player has won."""
        # Check horizontal
        for r in range(ROWS):
            for c in range(COLS - 3):
                if all(b[r][c+i] == player for i in range(4)):
                    return True
        # Check vertical
        for r in range(ROWS - 3):
            for c in range(COLS):
                if all(b[r+i][c] == player for i in range(4)):
                    return True
        # Check diagonal (down-right)
        for r in range(ROWS - 3):
            for c in range(COLS - 3):
                if all(b[r+i][c+i] == player for i in range(4)):
                    return True
        # Check diagonal (up-right)
        for r in range(3, ROWS):
            for c in range(COLS - 3):
                if all(b[r-i][c+i] == player for i in range(4)):
                    return True
        return False

    optimal_moves = []
    
    # Iterate through all possible first moves for 'O'
    for c in range(COLS):
        r = get_next_open_row(board, c)
        if r is None:
            continue

        board_after_my_move = copy.deepcopy(board)
        board_after_my_move[r][c] = ME

        # Check for an immediate win (win in 1)
        if check_win(board_after_my_move, ME):
            # This would be the best move, but let's assume none exist as per puzzle type
            optimal_moves.append((r, c))
            continue
            
        # Check for a forced win (win in 2)
        is_forced_win = True
        
        # Check every possible response from the opponent
        has_opponent_move = False
        for op_c in range(COLS):
            op_r = get_next_open_row(board_after_my_move, op_c)
            if op_r is None:
                continue

            has_opponent_move = True
            board_after_op_move = copy.deepcopy(board_after_my_move)
            board_after_op_move[op_r][op_c] = OPPONENT

            # After opponent's move, can I make a winning move?
            i_can_win = False
            for my_next_c in range(COLS):
                my_next_r = get_next_open_row(board_after_op_move, my_next_c)
                if my_next_r is None:
                    continue

                final_board = copy.deepcopy(board_after_op_move)
                final_board[my_next_r][my_next_c] = ME
                if check_win(final_board, ME):
                    i_can_win = True
                    break # Found a winning response
            
            # If there's even one opponent move for which I have no winning response,
            # my initial move was not a forced win.
            if not i_can_win:
                is_forced_win = False
                break
        
        # If the opponent has no moves and I haven't won, it's a draw/stalemate.
        if not has_opponent_move:
             is_forced_win = False

        if is_forced_win:
            optimal_moves.append((r, c))

    # Format the output
    col_map = "abcdefg"
    # The puzzle rows are 1-6 from top to bottom. Our board rows are 0-5.
    # Puzzle row = 6 - board_row_index
    formatted_moves = [f"{col_map[c]}{6 - r}" for r, c in optimal_moves]
    print(", ".join(sorted(formatted_moves)))

solve_connect_four()