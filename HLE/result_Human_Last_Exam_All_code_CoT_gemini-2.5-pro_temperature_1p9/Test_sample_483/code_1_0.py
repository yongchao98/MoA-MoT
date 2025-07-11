import sys

def solve_game():
    """
    Encapsulates and runs the game solver to find the shortest forced win for Player 1.
    """
    # It's a game tree search, a higher recursion limit might be needed for complex games.
    sys.setrecursionlimit(2000)

    # --- Game Configuration ---
    # Piece indices for the state tuple for easy access
    K1, N1, R1, R2, N2, K2 = 0, 1, 2, 3, 4, 5
    PLAYER1_PIECES = {K1, N1, R1}
    PLAYER2_PIECES = {K2, N2, R2}
    
    INITIAL_STATE = (0, 1, 2, 5, 6, 7)  # Positions: (K1, N1, R1, R2, N2, K2)
    BOARD_SIZE = 8
    CAPTURED = -1
    INFINITY = float('inf')

    # Memoization table for dynamic programming (state, player) -> result
    memo = {}

    # --- Helper Functions ---
    def get_board_map(state):
        """Creates a map from position -> piece_index for all occupied squares."""
        board_map = {}
        for piece_idx, pos in enumerate(state):
            if pos != CAPTURED:
                board_map[pos] = piece_idx
        return board_map

    def is_king_attacked_by_rook(state, king_idx, opponent_rook_idx):
        """Checks if a king is under direct attack by the opponent's rook."""
        king_pos = state[king_idx]
        rook_pos = state[opponent_rook_idx]

        if king_pos == CAPTURED or rook_pos == CAPTURED:
            return False

        board_map = get_board_map(state)
        start, end = min(king_pos, rook_pos), max(king_pos, rook_pos)

        # Check for any pieces between the king and the rook
        for i in range(start + 1, end):
            if i in board_map:
                return False  # Path is blocked, king is safe
        
        return True  # Path is clear, king is attacked

    def generate_legal_next_states(state, player):
        """Generates all possible legal next states for the current player."""
        legal_states = []
        my_pieces = PLAYER1_PIECES if player == 1 else PLAYER2_PIECES
        my_king_idx = K1 if player == 1 else K2
        opponent_rook_idx = R2 if player == 1 else R1

        board_map = get_board_map(state)
        my_positions = {state[p_idx] for p_idx in my_pieces if state[p_idx] != CAPTURED}

        for piece_idx in my_pieces:
            from_pos = state[piece_idx]
            if from_pos == CAPTURED:
                continue

            potential_destinations = []
            # King (K): Can move one step
            if piece_idx == K1 or piece_idx == K2:
                potential_destinations.extend([from_pos - 1, from_pos + 1])
            # Knight (N): Can move two steps
            elif piece_idx == N1 or piece_idx == N2:
                potential_destinations.extend([from_pos - 2, from_pos + 2])
            # Rook (R): Can move until blocked
            elif piece_idx == R1 or piece_idx == R2:
                # Move right
                for p in range(from_pos + 1, BOARD_SIZE):
                    potential_destinations.append(p)
                    if p in board_map: break
                # Move left
                for p in range(from_pos - 1, -1, -1):
                    potential_destinations.append(p)
                    if p in board_map: break

            for to_pos in potential_destinations:
                # Rule 1 & 2: Check if destination is valid and on board
                if not (0 <= to_pos < BOARD_SIZE): continue
                # Rule 3: Check for friendly fire
                if to_pos in my_positions: continue

                # Create a potential next state
                next_state_list = list(state)
                next_state_list[piece_idx] = to_pos
                
                # Handle capture
                if to_pos in board_map:
                    captured_piece_idx = board_map[to_pos]
                    next_state_list[captured_piece_idx] = CAPTURED
                
                next_state = tuple(next_state_list)

                # Rule 4: Final check for King's safety
                if not is_king_attacked_by_rook(next_state, my_king_idx, opponent_rook_idx):
                    legal_states.append(next_state)

        return legal_states

    # --- Minimax Solver ---
    def find_shortest_win_path(state, player):
        """
        Calculates the outcome of the game from the current state using minimax.
        """
        state_key = (state, player)
        if state_key in memo:
            return memo[state_key]

        # --- Base Cases ---
        if state[K2] == CAPTURED: return 0  # Player 1 wins
        if state[K1] == CAPTURED: return INFINITY # Player 1 loses

        legal_next_states = generate_legal_next_states(state, player)

        # Stalemate condition
        if not legal_next_states:
            return INFINITY

        # --- Recursive Step ---
        outcomes = [find_shortest_win_path(next_s, 3 - player) for next_s in legal_next_states]

        if player == 1:
            # P1 wants to win in the minimum number of moves
            result = 1 + min(outcomes)
        else: # player == 2
            # P2 wants to delay defeat for the maximum number of moves
            result = 1 + max(outcomes)

        memo[state_key] = result
        return result

    # --- Main Execution ---
    result = find_shortest_win_path(INITIAL_STATE, 1)
    print(int(result))

solve_game()