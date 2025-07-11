import math
from functools import lru_cache

def solve_game():
    """
    This function sets up and solves the game, printing the final answer.
    All helper functions and logic are defined within this scope.
    """

    # --- Game State and Rules ---

    # Piece indices in the state tuple
    K1, N1, R1, K2, N2, R2 = 0, 1, 2, 3, 4, 5
    P1_PIECES = (K1, N1, R1)
    P2_PIECES = (K2, N2, R2)

    # Player representation
    P1 = 1
    P2 = 2

    # Initial board configuration represented by piece positions
    # (pos_K1, pos_N1, pos_R1, pos_K2, pos_N2, pos_R2)
    INITIAL_STATE = (0, 1, 2, 7, 6, 5)

    # Outcome constants
    P1_WIN = 1
    P2_WIN = -1
    DRAW = 0

    # --- Helper and Move Generation Functions ---

    def get_board_map(board_tuple):
        """Creates a mapping from position to piece index for fast lookups."""
        board_map = {}
        for piece_idx, pos in enumerate(board_tuple):
            if pos != -1:  # -1 means captured
                board_map[pos] = piece_idx
        return board_map

    def is_king_in_check(board_tuple, player):
        """
        Checks if the specified player's king is under attack by the opponent's rook.
        A king is in check if there is an empty line of sight to the opponent's rook.
        """
        if player == P1:
            king_pos = board_tuple[K1]
            rook_pos = board_tuple[R2]
        else:  # player == P2
            king_pos = board_tuple[K2]
            rook_pos = board_tuple[R1]

        if king_pos == -1 or rook_pos == -1:
            return False

        board_map = get_board_map(board_tuple)
        start, end = sorted((king_pos, rook_pos))
        
        for i in range(start + 1, end):
            if i in board_map:
                return False  # Path is blocked
        
        return True

    def generate_legal_moves(board_tuple, player):
        """Generates all legal moves for a given player from a board state."""
        legal_moves = []
        board_map = get_board_map(board_tuple)
        
        my_pieces = P1_PIECES if player == P1 else P2_PIECES
        opponent_pieces = P2_PIECES if player == P1 else P1_PIECES
        
        for piece_idx in my_pieces:
            pos = board_tuple[piece_idx]
            if pos == -1:
                continue

            # Define potential moves based on piece type
            if piece_idx in (K1, K2):  # King
                potential_moves = [pos - 1, pos + 1]
            elif piece_idx in (N1, N2):  # Knight
                potential_moves = [pos - 2, pos + 2]
            elif piece_idx in (R1, R2):  # Rook
                potential_moves = []
                # Scan left
                for new_pos in range(pos - 1, -1, -1):
                    potential_moves.append(new_pos)
                    if new_pos in board_map: break
                # Scan right
                for new_pos in range(pos + 1, 8):
                    potential_moves.append(new_pos)
                    if new_pos in board_map: break
            
            for new_pos in potential_moves:
                if not (0 <= new_pos <= 7):
                    continue

                # A move is valid if destination is empty or occupied by an opponent
                if new_pos not in board_map or board_map[new_pos] in opponent_pieces:
                    new_board_list = list(board_tuple)
                    new_board_list[piece_idx] = new_pos
                    
                    # Handle capture
                    if new_pos in board_map:
                        captured_piece_idx = board_map[new_pos]
                        new_board_list[captured_piece_idx] = -1
                    
                    new_board_tuple = tuple(new_board_list)
                    
                    # A move is legal only if it does not leave one's own king in check
                    if not is_king_in_check(new_board_tuple, player):
                        legal_moves.append(new_board_tuple)
                        
        return legal_moves

    # --- Minimax Solver ---

    @lru_cache(maxsize=None)
    def solve(board_tuple, player):
        """
        Recursively solves the game using minimax with memoization.
        Returns: (outcome, plies_to_outcome)
        """
        moves = generate_legal_moves(board_tuple, player)

        if not moves:
            if is_king_in_check(board_tuple, player):
                # Checkmate: current player loses
                return (P2_WIN if player == P1 else P1_WIN), 0
            else:
                # Stalemate: draw
                return DRAW, 0

        child_results = [solve(move, P2 if player == P1 else P1) for move in moves]

        if player == P1:  # Player 1 is the maximizer
            wins = [res for res in child_results if res[0] == P1_WIN]
            if wins:
                min_plies = min(p for _, p in wins)
                return P1_WIN, 1 + min_plies
            
            if any(res[0] == DRAW for res in child_results):
                return DRAW, 100 # Plies to draw is high to prefer winning
            
            losses = [res for res in child_results if res[0] == P2_WIN]
            max_plies = max(p for _, p in losses)
            return P2_WIN, 1 + max_plies
        else:  # Player 2 is the minimizer
            wins = [res for res in child_results if res[0] == P2_WIN]
            if wins:
                min_plies = min(p for _, p in wins)
                return P2_WIN, 1 + min_plies
            
            if any(res[0] == DRAW for res in child_results):
                return DRAW, 100

            losses = [res for res in child_results if res[0] == P1_WIN]
            max_plies = max(p for _, p in losses)
            return P1_WIN, 1 + max_plies

    # --- Main Execution ---
    
    outcome, plies = solve(INITIAL_STATE, P1)
    
    if outcome == P1_WIN:
        # A "turn" involves a move by P1 and potentially a move by P2.
        # plies = 1 (P1 moves, wins) -> 1 turn
        # plies = 3 (P1, P2, P1 wins) -> 2 turns
        # Formula: math.ceil(plies / 2)
        turns = math.ceil(plies / 2)
        print(turns)
    elif outcome == P2_WIN:
        plies_to_loss = plies
        print(f"Player 1 cannot force a win. Player 2 wins in {math.ceil(plies_to_loss / 2)} turns.")
    else:
        print("The game ends in a draw with optimal play.")

solve_game()