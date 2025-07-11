import sys

# Set a higher recursion limit for deep game tree exploration.
sys.setrecursionlimit(2000)

# --- Global Definitions ---

# Piece identifiers for clarity and consistency.
K1, N1, R1 = 'K1', 'N1', 'R1'  # Player 1
K2, N2, R2 = 'K2', 'N2', 'R2'  # Player 2
EMPTY = '  '

P1_PIECES = {K1, N1, R1}
P2_PIECES = {K2, N2, R2}

# Memoization table to store results of solved game states.
memo = {}

def get_piece_positions(board):
    """Returns a dictionary mapping piece names to their board positions."""
    positions = {}
    for i, piece in enumerate(board):
        if piece != EMPTY:
            positions[piece] = i
    return positions

def is_king_attacked(board, player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    positions = get_piece_positions(board)
    
    king_piece = K1 if player == 1 else K2
    rook_piece = R2 if player == 1 else R1

    # If the king or the attacking rook is not on the board, no attack is possible.
    if king_piece not in positions or rook_piece not in positions:
        return False

    king_pos = positions[king_piece]
    rook_pos = positions[rook_piece]

    start, end = min(king_pos, rook_pos), max(king_pos, rook_pos)
    
    # Check for any blocking pieces between the king and rook.
    for i in range(start + 1, end):
        if board[i] != EMPTY:
            return False  # Path is blocked.
            
    return True  # Path is clear, the king is under attack.

def generate_moves(board, player):
    """Generates all legal next board states for the given player."""
    moves = []
    positions = get_piece_positions(board)
    player_pieces = P1_PIECES if player == 1 else P2_PIECES
    opponent_pieces = P2_PIECES if player == 1 else P1_PIECES

    for piece, pos in positions.items():
        if piece not in player_pieces:
            continue

        potential_destinations = []
        # King moves: one step left or right.
        if piece.startswith('K'):
            potential_destinations = [pos - 1, pos + 1]
        # Knight moves: two steps left or right.
        elif piece.startswith('N'):
            potential_destinations = [pos - 2, pos + 2]
        # Rook moves: any number of steps until blocked.
        elif piece.startswith('R'):
            # Move right
            for new_pos in range(pos + 1, 8):
                potential_destinations.append(new_pos)
                if board[new_pos] != EMPTY: break
            # Move left
            for new_pos in range(pos - 1, -1, -1):
                potential_destinations.append(new_pos)
                if board[new_pos] != EMPTY: break
        
        for new_pos in potential_destinations:
            # Check if the move is on the board and the destination is valid.
            if 0 <= new_pos < 8 and (board[new_pos] == EMPTY or board[new_pos] in opponent_pieces):
                new_board_list = list(board)
                new_board_list[pos] = EMPTY
                new_board_list[new_pos] = piece
                new_board = tuple(new_board_list)
                # A move is only legal if it does not leave one's own king in check.
                if not is_king_attacked(new_board, player):
                    moves.append(new_board)
                    
    return moves

def solve(board, player):
    """
    Recursively solves the game state using minimax with memoization.
    Returns (outcome, plies):
    - outcome: 1 for P1 win, 2 for P2 win, 0 for draw.
    - plies: Number of half-moves until the outcome.
    """
    state = (board, player)
    if state in memo:
        return memo[state]

    # Check for a win by capturing the opponent's king.
    opponent_king = K2 if player == 1 else K1
    if opponent_king not in board:
        # The previous player captured the king and won.
        return (3 - player, 0)

    legal_moves = generate_moves(board, player)

    # Base case: No legal moves available.
    if not legal_moves:
        if is_king_attacked(board, player):
            # Checkmate: current player loses.
            result = (3 - player, 0)
        else:
            # Stalemate: draw.
            result = (0, 0)
        memo[state] = result
        return result

    # --- Recursive Step (Minimax Logic) ---
    
    # Player 1 (maximizer) seeks a win in the fewest moves.
    if player == 1:
        best_outcome = -1  # Represents a loss for P1.
        min_win_plies, max_draw_plies, max_loss_plies = float('inf'), -1, -1
        
        for move in legal_moves:
            outcome, plies = solve(move, 2)
            if outcome == 1: # P1 win
                best_outcome = max(best_outcome, 1)
                min_win_plies = min(min_win_plies, plies)
            elif outcome == 0: # Draw
                best_outcome = max(best_outcome, 0)
                max_draw_plies = max(max_draw_plies, plies)
            else: # P1 loss
                max_loss_plies = max(max_loss_plies, plies)

        if best_outcome == 1:
            result = (1, min_win_plies + 1)
        elif best_outcome == 0:
            result = (0, max_draw_plies + 1)
        else:
            result = (2, max_loss_plies + 1)
            
    # Player 2 (minimizer) seeks to avoid loss, or prolong it if inevitable.
    else: # player == 2
        best_outcome = 3 # Represents a loss for P2.
        min_win_plies, max_draw_plies, max_loss_plies = float('inf'), -1, -1
        
        for move in legal_moves:
            outcome, plies = solve(move, 1)
            if outcome == 2: # P2 win
                best_outcome = min(best_outcome, 2)
                min_win_plies = min(min_win_plies, plies)
            elif outcome == 0: # Draw
                best_outcome = min(best_outcome, 0)
                max_draw_plies = max(max_draw_plies, plies)
            else: # P2 loss
                max_loss_plies = max(max_loss_plies, plies)

        if best_outcome == 2:
            result = (2, min_win_plies + 1)
        elif best_outcome == 0:
            result = (0, max_draw_plies + 1)
        else: # All moves lead to a loss for P2, so stall.
            result = (1, max_loss_plies + 1)

    memo[state] = result
    return result

def main():
    """Sets up the initial board and starts the solving process."""
    initial_board = (K1, N1, R1, EMPTY, EMPTY, R2, N2, K2)
    
    # Start the game from Player 1's perspective.
    outcome, plies = solve(initial_board, 1)
    
    if outcome == 1:
        # A turn is one move by Player 1.
        # Plies 1 -> Turn 1
        # Plies 3 -> Turn 2
        # Plies 5 -> Turn 3
        # turns = (plies + 1) / 2
        turns = (plies + 1) // 2
        print(f"{turns}")
    elif outcome == 2:
        print("Player 2 can force a win.")
    else:
        print("The game is a draw with optimal play.")

if __name__ == '__main__':
    main()
    print("<<<4>>>")