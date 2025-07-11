import functools

# Use lru_cache for memoization to store results of solved states
@functools.lru_cache(maxsize=None)
def solve(board, player_to_move):
    """
    Recursively solves the game state using minimax.
    Returns a tuple: (status, plies).
    - status: 1 if the current player can force a win,
             -1 if the current player is forced to lose,
              0 if the game is a forced draw.
    - plies: The number of half-moves (plies) until the game ends.
    """
    # --- Base Cases: Check for game termination ---
    p1_king_pos = -1
    p2_king_pos = -1
    for i, p in enumerate(board):
        if p == 'K1': p1_king_pos = i
        elif p == 'K2': p2_king_pos = i

    if p2_king_pos == -1:  # P2's King is captured
        return (1, 0) if player_to_move == 'P1' else (-1, 0)
    if p1_king_pos == -1:  # P1's King is captured
        return (-1, 0) if player_to_move == 'P1' else (1, 0)

    legal_moves = generate_legal_moves(board, player_to_move)

    if not legal_moves:
        if is_king_in_check(board, player_to_move):
            # Checkmate: current player loses
            return -1, 0
        else:
            # Stalemate: draw
            return 0, 0

    # --- Recursive Step: Explore child states ---
    
    best_outcome = -1  # Assume loss
    best_plies = -1

    # These will store the results for different outcomes
    win_results = []
    draw_results = []
    loss_results = []

    next_player = 'P2' if player_to_move == 'P1' else 'P1'
    for move in legal_moves:
        status, plies = solve(move, next_player)
        # The status is from the opponent's perspective, so we flip it
        outcome = -status
        
        if outcome == 1:
            win_results.append(plies)
        elif outcome == 0:
            draw_results.append(plies)
        else: # outcome == -1
            loss_results.append(plies)

    if win_results:
        # If there's a path to a win, choose the one with the fewest plies.
        best_outcome = 1
        best_plies = min(win_results)
    elif draw_results:
        # If no win, but can force a draw, take it.
        best_outcome = 0
        # In a draw, we don't care about plies, but let's be consistent.
        best_plies = max(draw_results) if draw_results else 0
    else:
        # If all moves lead to a loss, choose the one that stalls the longest.
        best_outcome = -1
        best_plies = max(loss_results)

    return best_outcome, best_plies + 1

def is_king_in_check(board, player):
    """Checks if the specified player's King is under attack."""
    king_piece = 'K1' if player == 'P1' else 'K2'
    rook_piece = 'R2' if player == 'P1' else 'R1'

    king_pos = -1
    rook_pos = -1
    for i, p in enumerate(board):
        if p == king_piece: king_pos = i
        elif p == rook_piece: rook_pos = i
    
    if king_pos == -1 or rook_pos == -1:
        return False

    # Check for clear line of sight
    start, end = sorted((king_pos, rook_pos))
    for i in range(start + 1, end):
        if board[i] is not None:
            return False # Path is blocked
    return True

def generate_legal_moves(board, player):
    """Generates all legal moves for the given player."""
    legal_moves = []
    player_pieces = ('K1', 'N1', 'R1') if player == 'P1' else ('K2', 'N2', 'R2')

    for pos, piece in enumerate(board):
        if piece not in player_pieces:
            continue

        potential_moves = []
        # --- King Moves ---
        if piece.startswith('K'):
            for d in [-1, 1]:
                new_pos = pos + d
                if 0 <= new_pos <= 7 and (board[new_pos] is None or board[new_pos] not in player_pieces):
                    potential_moves.append(new_pos)
        
        # --- Knight Moves ---
        elif piece.startswith('N'):
            for d in [-2, 2]:
                new_pos = pos + d
                if 0 <= new_pos <= 7 and (board[new_pos] is None or board[new_pos] not in player_pieces):
                    potential_moves.append(new_pos)

        # --- Rook Moves ---
        elif piece.startswith('R'):
            for d in [-1, 1]: # Direction
                for i in range(1, 8):
                    new_pos = pos + d * i
                    if not (0 <= new_pos <= 7):
                        break # Off board
                    if board[new_pos] is not None:
                        if board[new_pos] in player_pieces:
                            break # Blocked by own piece
                        else:
                            potential_moves.append(new_pos) # Capture
                            break # Cannot move further
                    potential_moves.append(new_pos) # Empty square

        # Validate each potential move for king safety
        for new_pos in potential_moves:
            new_board = list(board)
            new_board[new_pos] = piece
            new_board[pos] = None
            new_board_tuple = tuple(new_board)
            if not is_king_in_check(new_board_tuple, player):
                legal_moves.append(new_board_tuple)
                
    return legal_moves

def main():
    """
    Sets up the initial game state and calls the solver.
    """
    initial_board = (
        'K1', 'N1', 'R1', None, None, 'R2', 'N2', 'K2'
    )
    
    # Solve from the perspective of Player 1
    status, plies = solve(initial_board, 'P1')

    if status == 1:
        # A turn consists of a move by P1.
        # Plies are half-moves.
        # 1 ply -> 1 turn
        # 3 plies -> 2 turns
        # 5 plies -> 3 turns
        # 7 plies -> 4 turns
        # Formula: ceil(plies / 2) or (plies + 1) // 2
        p1_moves = (plies + 1) // 2
        print(f"{p1_moves}")
    elif status == 0:
        print("The game is a forced draw.")
    else:
        print("Player 1 cannot force a win.")

if __name__ == "__main__":
    main()