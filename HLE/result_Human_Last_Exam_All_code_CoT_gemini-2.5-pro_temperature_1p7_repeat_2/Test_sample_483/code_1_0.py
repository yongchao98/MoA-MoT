import sys

# It's recommended to increase the recursion limit for deep game tree searches.
sys.setrecursionlimit(3000)

# Global cache for memoization
memo = {}
# For cycle detection in the recursion path
recursion_path = set()

# Piece indices for the positions tuple
K1, N1, R1, K2, N2, R2 = 0, 1, 2, 3, 4, 5
PIECE_NAMES = {K1: "K1", N1: "N1", R1: "R1", K2: "K2", N2: "N2", R2: "R2"}


# Board boundaries
BOARD_MIN, BOARD_MAX = 0, 7

def get_board_map(positions):
    """Creates a map from position to piece index for quick lookups."""
    board_map = {}
    for piece_idx, pos in enumerate(positions):
        if pos is not None:
            board_map[pos] = piece_idx
    return board_map

def is_king_attacked_by_rook(positions, player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    if player == 1:
        king_idx, rook_idx = K1, R2
    else: # player == 2
        king_idx, rook_idx = K2, R1

    king_pos = positions[king_idx]
    rook_pos = positions[rook_idx]

    if king_pos is None or rook_pos is None:
        return False

    board_map = get_board_map(positions)
    start = min(king_pos, rook_pos) + 1
    end = max(king_pos, rook_pos)
    
    for pos in range(start, end):
        if pos in board_map:
            return False  # Path is blocked
    
    return True

def get_legal_moves(positions, player):
    """Generates all legal moves for a given player."""
    legal_moves = []
    
    if player == 1:
        player_piece_indices = [K1, N1, R1]
    else: # player == 2
        player_piece_indices = [K2, N2, R2]

    board_map = get_board_map(positions)

    for piece_idx in player_piece_indices:
        current_pos = positions[piece_idx]
        if current_pos is None:
            continue
        
        # --- Generate potential destinations based on piece type ---
        potential_dests = []
        if piece_idx in [K1, K2]: # King
            potential_dests.extend([current_pos - 1, current_pos + 1])
        elif piece_idx in [N1, N2]: # Knight
            potential_dests.extend([current_pos - 2, current_pos + 2])
        elif piece_idx in [R1, R2]: # Rook
            # Move left until blocked or edge
            for dest in range(current_pos - 1, BOARD_MIN - 1, -1):
                potential_dests.append(dest)
                if dest in board_map: break
            # Move right until blocked or edge
            for dest in range(current_pos + 1, BOARD_MAX + 1):
                potential_dests.append(dest)
                if dest in board_map: break

        # --- Validate each potential destination ---
        for dest in potential_dests:
            # Check board boundaries
            if not (BOARD_MIN <= dest <= BOARD_MAX):
                continue
            
            # Check if destination is occupied by a friendly piece
            if dest in board_map and board_map[dest] in player_piece_indices:
                continue

            # Create a temporary new state to check for king safety
            next_positions = list(positions)
            next_positions[piece_idx] = dest
            
            # Handle capture for the check
            if dest in board_map:
                captured_piece_idx = board_map[dest]
                next_positions[captured_piece_idx] = None
            
            # A move is illegal if it places the player's own king in check
            if not is_king_attacked_by_rook(tuple(next_positions), player):
                legal_moves.append((piece_idx, dest))

    return legal_moves

def solve(positions, player):
    """
    Recursively solves the game state using minimax.
    Returns: (winner, turns_to_win), where turns_to_win is the number of plies.
    """
    state_key = (positions, player)
    
    # --- Base Cases, Memoization, and Cycle Detection ---
    if positions[K2] is None: return (1, 0) # P1 wins
    if positions[K1] is None: return (2, 0) # P2 wins
    if state_key in memo: return memo[state_key]
    if state_key in recursion_path: return (0, 0) # Draw by repetition

    recursion_path.add(state_key)

    moves = get_legal_moves(positions, player)
    
    if not moves:
        recursion_path.remove(state_key)
        # No legal moves means stalemate, which is a draw.
        memo[state_key] = (0, 0)
        return (0, 0)

    # --- Recursive Step ---
    outcomes = []
    for piece_idx, dest in moves:
        next_positions_list = list(positions)
        board_map = get_board_map(positions)
        
        # Apply move and handle capture
        next_positions_list[piece_idx] = dest
        if dest in board_map:
            captured_idx = board_map[dest]
            next_positions_list[captured_idx] = None
            
        outcomes.append(solve(tuple(next_positions_list), 3 - player))

    # --- Minimax Logic: determine the best outcome ---
    my_result = None
    if player == 1:
        # P1 looks for a win in min turns
        p1_wins = [r for r in outcomes if r[0] == 1]
        if p1_wins:
            min_turns = min(r[1] for r in p1_wins)
            my_result = (1, min_turns + 1)
        else:
            # If no win, prefer a draw
            if any(r[0] == 0 for r in outcomes):
                my_result = (0, 0)
            # If all moves lead to a loss, stall as long as possible
            else:
                p2_losses = [r for r in outcomes if r[0] == 2]
                max_turns = max(r[1] for r in p2_losses)
                my_result = (2, max_turns + 1)
    else: # player == 2
        # P2 looks for a win in min turns
        p2_wins = [r for r in outcomes if r[0] == 2]
        if p2_wins:
            min_turns = min(r[1] for r in p2_wins)
            my_result = (2, min_turns + 1)
        else:
            # If no win, prefer a draw
            if any(r[0] == 0 for r in outcomes):
                my_result = (0, 0)
            # If all moves lead to a loss, stall as long as possible
            else:
                p1_losses = [r for r in outcomes if r[0] == 1]
                max_turns = max(r[1] for r in p1_losses)
                my_result = (1, max_turns + 1)

    recursion_path.remove(state_key)
    memo[state_key] = my_result
    return my_result

if __name__ == "__main__":
    # Initial Configuration: [K1][N1][R1][ ][ ][R2][N2][K2]
    # Player 1: K1 at 0, N1 at 1, R1 at 2
    # Player 2: R2 at 5, N2 at 6, K2 at 7
    # The tuple represents positions of (K1, N1, R1, K2, N2, R2)
    initial_positions = (0, 1, 2, 7, 6, 5) 
    
    winner, turns = solve(initial_positions, 1)

    if winner == 1:
        print(f"Player 1 can force a win in {turns} turns (plies).")
        print(turns)
    elif winner == 2:
        print(f"Player 2 can force a win in {turns} turns (plies).")
        print(turns)
    else:
        print("The game results in a draw with optimal play.")
