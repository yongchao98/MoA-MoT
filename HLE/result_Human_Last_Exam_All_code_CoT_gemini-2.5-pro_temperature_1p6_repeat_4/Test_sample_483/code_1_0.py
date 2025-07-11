import math

# Using a global dictionary for memoization (transposition table)
memo = {}
PIECE_ORDER = ['K1', 'N1', 'R1', 'K2', 'N2', 'R2']
P1_PIECES = {'K1', 'N1', 'R1'}
P2_PIECES = {'K2', 'N2', 'R2'}

def get_board_and_pos_map(state_tuple):
    """Helper to convert state tuple to a board list and a piece->position map."""
    board = [None] * 8
    pos_map = {}
    for i, pos in enumerate(state_tuple):
        if pos is not None:
            piece_name = PIECE_ORDER[i]
            board[pos] = piece_name
            pos_map[piece_name] = pos
    return board, pos_map

def is_in_check(state_tuple, king_owner_player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    board, pos_map = get_board_and_pos_map(state_tuple)

    if king_owner_player == 'P1':
        king_pos = pos_map.get('K1')
        rook_pos = pos_map.get('R2')
    else:  # 'P2'
        king_pos = pos_map.get('K2')
        rook_pos = pos_map.get('R1')

    # If king or opponent's rook is captured, not in check.
    if king_pos is None or rook_pos is None:
        return False

    # Check for any pieces between the king and the rook.
    start, end = sorted((king_pos, rook_pos))
    for i in range(start + 1, end):
        if board[i] is not None:
            return False  # Path is blocked

    return True  # Path is clear, king is in check.

def generate_legal_moves(state_tuple, player):
    """Generates all legal next state tuples for the given player."""
    legal_next_states = []
    board, pos_map = get_board_and_pos_map(state_tuple)
    
    player_pieces = P1_PIECES if player == 'P1' else P2_PIECES
    king_to_protect = 'K1' if player == 'P1' else 'K2'

    for piece, start_pos in pos_map.items():
        if piece not in player_pieces:
            continue

        piece_type = piece[0]
        potential_end_positions = []

        if piece_type == 'K':  # King
            potential_end_positions.extend([start_pos - 1, start_pos + 1])
        elif piece_type == 'N':  # Knight
            potential_end_positions.extend([start_pos - 2, start_pos + 2])
        elif piece_type == 'R':  # Rook
            # Scan right
            for i in range(start_pos + 1, 8):
                potential_end_positions.append(i)
                if board[i] is not None: break
            # Scan left
            for i in range(start_pos - 1, -1, -1):
                potential_end_positions.append(i)
                if board[i] is not None: break
        
        for end_pos in potential_end_positions:
            if not (0 <= end_pos < 8):
                continue
            
            target_piece = board[end_pos]
            if target_piece and target_piece in player_pieces:
                continue

            # Create the potential next state
            next_state_list = list(state_tuple)
            next_state_list[PIECE_ORDER.index(piece)] = end_pos
            if target_piece:
                next_state_list[PIECE_ORDER.index(target_piece)] = None
            
            next_state_tuple = tuple(next_state_list)

            # A move is legal only if it does not leave the player's own king in check
            if not is_in_check(next_state_tuple, player):
                legal_next_states.append(next_state_tuple)
                
    return legal_next_states

def solve_game(state_tuple, player):
    """
    Recursively solves the game using minimax, returning the number of ply to a P1 win.
    P1 (minimizer) seeks the smallest ply count.
    P2 (maximizer) seeks the largest ply count.
    Returns float('inf') for a loss or draw for P1.
    """
    state_key = (state_tuple, player)
    if state_key in memo:
        return memo[state_key]

    # Check for terminal states (win/loss by king capture)
    if state_tuple[PIECE_ORDER.index('K2')] is None: return 0  # P1 wins
    if state_tuple[PIECE_ORDER.index('K1')] is None: return float('inf') # P2 wins

    # Generate moves and check for stalemate/checkmate
    legal_moves = generate_legal_moves(state_tuple, player)
    if not legal_moves:
        if is_in_check(state_tuple, player):
            return float('inf') # Checkmate (current player loses)
        else:
            return float('inf') # Stalemate (draw)

    if player == 'P1':  # P1 wants to MINIMIZE the ply to a win
        min_ply = float('inf')
        for next_state in legal_moves:
            ply = solve_game(next_state, 'P2')
            min_ply = min(min_ply, ply)
        result = 1 + min_ply
    else:  # 'P2', P2 wants to MAXIMIZE the ply to a P1 win
        max_ply = -float('inf')
        for next_state in legal_moves:
            ply = solve_game(next_state, 'P1')
            max_ply = max(max_ply, ply)
        result = 1 + max_ply if max_ply != float('inf') else float('inf')

    memo[state_key] = result
    return result

def main():
    """Main function to solve the game and print the result."""
    # Initial state: (K1, N1, R1, K2, N2, R2)
    initial_state = (0, 1, 2, 7, 6, 5)
    
    # Solve from P1's perspective
    total_ply = solve_game(initial_state, 'P1')

    if total_ply == float('inf'):
        print("Player 1 cannot force a win.")
    else:
        # Number of P1 turns is ceil(ply / 2)
        p1_turns = (int(total_ply) + 1) // 2
        print("Forced win analysis complete.")
        print(f"Total ply for forced win: {int(total_ply)}")
        print(f"Calculation for Player 1 turns: ({int(total_ply)} + 1) // 2 = {p1_turns}")
        print(f"Player 1 can force a win in {p1_turns} turns.")
        print(f"<<<{p1_turns}>>>")

if __name__ == "__main__":
    main()
