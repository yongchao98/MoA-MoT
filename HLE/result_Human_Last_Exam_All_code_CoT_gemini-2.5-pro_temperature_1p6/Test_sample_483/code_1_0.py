import sys

# Set a higher recursion limit for deep game trees, though it may not be necessary for this problem.
sys.setrecursionlimit(2000)

# Memoization table to cache results of states that have already been evaluated.
# Key: (board_tuple, player_num), Value: (outcome, ply)
memo = {}

# --- Helper Functions ---

def find_piece(board, piece_name):
    """Finds the position of a piece on the board."""
    try:
        return board.index(piece_name)
    except ValueError:
        return None

def is_king_in_check(board, player_num):
    """
    Checks if the specified player's King is under attack by the opponent's Rook.
    A King is in check if the Rook has a clear line of sight to it.
    """
    king_piece = f"K{player_num}"
    opponent_rook = f"R{2 if player_num == 1 else 1}"

    king_pos = find_piece(board, king_piece)
    rook_pos = find_piece(board, opponent_rook)

    # If either piece is not on the board, the king is not in check by that piece.
    if king_pos is None or rook_pos is None:
        return False

    start, end = min(king_pos, rook_pos), max(king_pos, rook_pos)
    # Check if the path between the King and Rook is empty.
    for i in range(start + 1, end):
        if board[i] != " ":
            return False  # Path is blocked.
    
    return True  # Path is clear, King is in check.

def get_legal_moves(board, player_num):
    """
    Generates all legal moves for a given player.
    A move is legal if it follows the piece's movement rules, the destination is
    valid, and it does not leave the player's own king in check.
    """
    legal_moves = []
    player_char = str(player_num)
    
    for pos, piece in enumerate(board):
        if player_char not in piece:
            continue

        targets = []
        piece_type = piece[0]

        # Generate potential target positions based on piece type.
        if piece_type == 'K':
            targets = [pos - 1, pos + 1]
        elif piece_type == 'N':
            targets = [pos - 2, pos + 2]
        elif piece_type == 'R':
            # Move left
            for i in range(pos - 1, -1, -1):
                targets.append(i)
                if board[i] != " ": break
            # Move right
            for i in range(pos + 1, 8):
                targets.append(i)
                if board[i] != " ": break

        for target_pos in targets:
            # Rule 1: Destination must be within the board.
            if not (0 <= target_pos < 8):
                continue
            
            # Rule 2: Destination cannot be occupied by a friendly piece.
            if player_char in board[target_pos]:
                continue
            
            # Simulate the move on a temporary board.
            new_board = list(board)
            new_board[target_pos] = piece
            new_board[pos] = " "
            
            # Rule 3: The move must not result in the player's own king being in check.
            if not is_king_in_check(new_board, player_num):
                legal_moves.append(new_board)
                
    return legal_moves

def solve(board_tuple, player_num):
    """
    Recursively solves the game from a given state using minimax.
    Returns:
        - outcome: 1 for P1 win, -1 for P2 win (P1 loss), 0 for draw.
        - ply: The number of half-moves until the game ends.
    """
    # Check memoization table first
    if (board_tuple, player_num) in memo:
        return memo[(board_tuple, player_num)]

    # --- Base Cases (Terminal States) ---

    # Win/Loss by King capture
    if find_piece(board_tuple, "K1") is None:
        return (-1, 0)  # Player 1's King is captured, P1 loses.
    if find_piece(board_tuple, "K2") is None:
        return (1, 0)   # Player 2's King is captured, P1 wins.

    # Generate all legal moves from the current position.
    legal_moves = get_legal_moves(list(board_tuple), player_num)

    # Stalemate or Checkmate
    if not legal_moves:
        if is_king_in_check(list(board_tuple), player_num):
            # No legal moves and king is in check: Checkmate. The current player loses.
            result = (-1, 0) if player_num == 1 else (1, 0)
        else:
            # No legal moves and king is not in check: Stalemate.
            result = (0, 0)
        memo[(board_tuple, player_num)] = result
        return result

    # --- Recursive Step (Minimax) ---
    
    # Recursively call solve for all possible next states.
    opponent_num = 2 if player_num == 1 else 1
    child_results = [solve(tuple(move), opponent_num) for move in legal_moves]

    if player_num == 1:  # Player 1 (Maximizer)
        wins = [(o, p) for o, p in child_results if o == 1]
        draws = [(o, p) for o, p in child_results if o == 0]
        
        if wins:
            # If there are winning moves, choose the one with the minimum ply (fastest win).
            best_outcome = 1
            best_ply = min(p for _, p in wins) + 1
        elif draws:
            # If no wins but can draw, choose the drawing line that lasts longest.
            best_outcome = 0
            best_ply = max(p for _, p in draws) + 1
        else:
            # If all moves lead to a loss, choose the one with maximum ply (stall).
            best_outcome = -1
            best_ply = max(p for _, p in child_results) + 1
    else:  # Player 2 (Minimizer)
        losses = [(o, p) for o, p in child_results if o == -1] # A P1 loss is a P2 win
        draws = [(o, p) for o, p in child_results if o == 0]
        
        if losses:
            # If P2 can win, choose the fastest win.
            best_outcome = -1
            best_ply = min(p for _, p in losses) + 1
        elif draws:
            # If no wins for P2 but can draw, choose longest draw.
            best_outcome = 0
            best_ply = max(p for _, p in draws) + 1
        else:
            # If all moves lead to P2 losing, stall as long as possible.
            best_outcome = 1
            best_ply = max(p for _, p in child_results) + 1
            
    # Cache and return the result for the current state.
    memo[(board_tuple, player_num)] = (best_outcome, best_ply)
    return best_outcome, best_ply

if __name__ == "__main__":
    # Initial board configuration
    initial_board = ("K1", "N1", "R1", " ", " ", "R2", "N2", "K2")
    
    # Start the analysis from the initial state, with Player 1 to move.
    outcome, ply = solve(initial_board, 1)
    
    # A turn consists of a move by P1 and a move by P2.
    # We use integer division to calculate turns from ply.
    turns = (ply + 1) // 2
    
    # Final output as requested
    print(f"Outcome: Player 1 {'wins' if outcome == 1 else ('loses' if outcome == -1 else 'draws')}")
    print(f"Number of ply (half-moves) to force result: {ply}")
    print("Calculation for turns:")
    print(f"({ply} + 1) // 2 = {turns}")
    print(f"\nPlayer 1 can force a win in {turns} turns.")
    print(f"<<<{turns}>>>")