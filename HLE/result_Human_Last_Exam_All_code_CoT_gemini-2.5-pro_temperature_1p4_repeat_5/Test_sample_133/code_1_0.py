import copy

def solve():
    """
    Solves the Connect 4 puzzle to find the optimal moves for 'O'.
    """
    board = [
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', 'O', '.', '.'],
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    ]
    player = 'O'
    opponent = 'X'
    rows, cols = 6, 7

    def get_lowest_empty_row(board, col):
        """Finds the lowest empty row in a given column."""
        for r in range(rows - 1, -1, -1):
            if board[r][col] == '.':
                return r
        return -1

    def get_valid_moves(board):
        """Gets all valid moves."""
        valid_moves = []
        for c in range(cols):
            r = get_lowest_empty_row(board, c)
            if r != -1:
                valid_moves.append((r, c))
        return valid_moves

    def check_win(board, player):
        """Checks if the specified player has won."""
        # Check horizontal
        for r in range(rows):
            for c in range(cols - 3):
                if all(board[r][c+i] == player for i in range(4)):
                    return True
        # Check vertical
        for r in range(rows - 3):
            for c in range(cols):
                if all(board[r+i][c] == player for i in range(4)):
                    return True
        # Check diagonal (down-right)
        for r in range(rows - 3):
            for c in range(cols - 3):
                if all(board[r+i][c+i] == player for i in range(4)):
                    return True
        # Check diagonal (up-right)
        for r in range(3, rows):
            for c in range(cols - 3):
                if all(board[r-i][c+i] == player for i in range(4)):
                    return True
        return False

    def find_winning_spots(board, player):
        """Finds all spots where the player can play to win immediately."""
        winning_spots = set()
        for r, c in get_valid_moves(board):
            temp_board = copy.deepcopy(board)
            temp_board[r][c] = player
            if check_win(temp_board, player):
                winning_spots.add((r, c))
        return winning_spots

    def to_connect4_notation(r, c):
        """Converts (row, col) to Connect 4 notation like 'a1'."""
        col_char = chr(ord('a') + c)
        row_char = str(r + 1)
        return f"{col_char}{row_char}"

    # --- Main Logic ---

    # Level 1: Find moves that win immediately.
    winning_moves_now = find_winning_spots(board, player)
    if winning_moves_now:
        print(", ".join(sorted([to_connect4_notation(r, c) for r, c in winning_moves_now])))
        return

    # Level 2: Find moves that create a fork (unblockable win in 2).
    fork_moves = []
    for r_move, c_move in get_valid_moves(board):
        # Make the move for 'O' on a temporary board
        board_after_O_move = copy.deepcopy(board)
        board_after_O_move[r_move][c_move] = player

        # Find where 'O' can win on the next turn
        o_winning_spots = find_winning_spots(board_after_O_move, player)

        # If there are fewer than 2 winning spots, it's not a fork.
        if len(o_winning_spots) < 2:
            continue

        # Check if opponent can block one threat, but another remains.
        can_be_blocked = False
        for r_block, c_block in o_winning_spots:
            # Simulate opponent blocking one of the threats
            board_after_X_block = copy.deepcopy(board_after_O_move)
            # The row to block might be different if the column isn't empty up to r_block
            block_row_actual = get_lowest_empty_row(board_after_O_move, c_block)
            if block_row_actual == -1: continue
            
            board_after_X_block[block_row_actual][c_block] = opponent

            # Check if 'O' still has a winning move
            if not find_winning_spots(board_after_X_block, player):
                can_be_blocked = True
                break
        
        if not can_be_blocked:
            fork_moves.append((r_move, c_move))
            
    if fork_moves:
        print(", ".join(sorted([to_connect4_notation(r, c) for r, c in fork_moves])))
        return

    print("No optimal winning moves found.")


solve()