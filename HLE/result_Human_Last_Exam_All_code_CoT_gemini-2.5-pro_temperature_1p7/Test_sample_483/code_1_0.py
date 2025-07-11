import functools

def solve_game():
    """
    This script solves the described board game to find the minimum number of turns
    for Player 1 to force a win against an optimal opponent.
    """

    # --- Game State and Rules Implementation ---

    # Use integers to represent pieces for efficiency
    K1, N1, R1, K2, N2, R2 = 0, 1, 2, 3, 4, 5
    PLAYER_1_PIECES = {K1, N1, R1}
    PLAYER_2_PIECES = {K2, N2, R2}
    PIECE_NAMES = ["K1", "N1", "R1", "K2", "N2", "R2"]
    BOARD_SIZE = 8

    def get_board_from_positions(positions):
        """Helper to create a list representation of the board."""
        board = [None] * BOARD_SIZE
        for piece_idx, pos in enumerate(positions):
            if pos is not None:
                board[pos] = piece_idx
        return board

    def is_king_in_check(positions, player):
        """Checks if the specified player's King is under attack by the opponent's Rook."""
        if player == 1:
            king_pos, rook_pos = positions[K1], positions[R2]
        else: # player == 2
            king_pos, rook_pos = positions[K2], positions[R1]

        if king_pos is None or rook_pos is None:
            return False

        board = get_board_from_positions(positions)
        start, end = sorted((king_pos, rook_pos))
        
        # Check for any pieces between the king and the rook
        for i in range(start + 1, end):
            if board[i] is not None:
                return False  # Path is blocked
        
        return True

    @functools.lru_cache(maxsize=None)
    def generate_legal_moves(state):
        """Generates all legal moves for the current player in a given state."""
        positions, turn = state
        legal_moves = []
        player_pieces = PLAYER_1_PIECES if turn == 1 else PLAYER_2_PIECES
        board = get_board_from_positions(positions)

        for piece_idx in player_pieces:
            pos = positions[piece_idx]
            if pos is None:
                continue
            
            piece_type = PIECE_NAMES[piece_idx][0]
            possible_dests = []

            if piece_type == 'K':
                if pos > 0: possible_dests.append(pos - 1)
                if pos < BOARD_SIZE - 1: possible_dests.append(pos + 1)
            elif piece_type == 'N':
                if pos > 1: possible_dests.append(pos - 2)
                if pos < BOARD_SIZE - 2: possible_dests.append(pos + 2)
            elif piece_type == 'R':
                # Move left
                for i in range(pos - 1, -1, -1):
                    possible_dests.append(i)
                    if board[i] is not None: break
                # Move right
                for i in range(pos + 1, BOARD_SIZE):
                    possible_dests.append(i)
                    if board[i] is not None: break
            
            for dest in possible_dests:
                # Rule: Cannot land on a friendly piece
                if board[dest] is not None and board[dest] in player_pieces:
                    continue

                # Simulate the move to check for king safety
                temp_positions = list(positions)
                if board[dest] is not None:
                    captured_piece_idx = board[dest]
                    temp_positions[captured_piece_idx] = None
                temp_positions[piece_idx] = dest
                
                # Rule: King cannot be in check after the move
                if not is_king_in_check(tuple(temp_positions), turn):
                    legal_moves.append((piece_idx, dest))
        return legal_moves

    # --- Minimax Solver with Memoization ---

    @functools.lru_cache(maxsize=None)
    def solve_state(state):
        """
        Recursively determines the outcome of the game from the given state.
        Returns a tuple: (outcome_string, plies_to_outcome)
        """
        positions, turn = state
        moves = generate_legal_moves(state)

        # Base Case: No legal moves available
        if not moves:
            if is_king_in_check(positions, turn): # Checkmate
                return ("P2_WIN", 0) if turn == 1 else ("P1_WIN", 0)
            else: # Stalemate
                return ("DRAW", float('inf'))

        outcomes = []
        for piece_idx, dest in moves:
            next_positions_list = list(positions)
            board = get_board_from_positions(positions)

            # Apply capture if destination is occupied by opponent
            if board[dest] is not None:
                captured_piece_idx = board[dest]
                next_positions_list[captured_piece_idx] = None
            
            next_positions_list[piece_idx] = dest
            next_positions = tuple(next_positions_list)

            # Check for immediate win by capturing the King
            if turn == 1 and next_positions[K2] is None:
                outcomes.append(("P1_WIN", 1))
                continue
            if turn == 2 and next_positions[K1] is None:
                outcomes.append(("P2_WIN", 1))
                continue

            # Recurse on the next state
            next_turn = 2 if turn == 1 else 1
            result, plies = solve_state((next_positions, next_turn))
            outcomes.append((result, plies + 1))

        # Apply Minimax logic based on the current player
        if turn == 1:  # Player 1 wants to win as fast as possible
            wins = [p for res, p in outcomes if res == "P1_WIN"]
            if wins: return ("P1_WIN", min(wins))
            draws = [p for res, p in outcomes if res == "DRAW"]
            if draws: return ("DRAW", float('inf'))
            losses = [p for res, p in outcomes if res == "P2_WIN"]
            return ("P2_WIN", max(losses)) # If losing, prolong the game
        else:  # Player 2 wants to win fast or stall a loss
            wins = [p for res, p in outcomes if res == "P2_WIN"]
            if wins: return ("P2_WIN", min(wins))
            draws = [p for res, p in outcomes if res == "DRAW"]
            if draws: return ("DRAW", float('inf'))
            losses = [p for res, p in outcomes if res == "P1_WIN"]
            return ("P1_WIN", max(losses)) # If losing, prolong the game

    # --- Main Execution ---
    
    print("Solving the game, please wait...")
    
    # Initial state: [K1][N1][R1][ ][ ][R2][N2][K2]
    # Positions are 0-indexed: (K1, N1, R1, K2, N2, R2)
    initial_positions = (0, 1, 2, 7, 6, 5)
    initial_state = (initial_positions, 1) # Player 1's turn

    result, plies = solve_state(initial_state)

    if result == "P1_WIN":
        # A "ply" is one player's move. A "turn" involves a move from both players.
        # Plies  Turns (for P1 win)
        # 1      1  (P1 wins on move 1)
        # 3      2  (P1 wins on move 2)
        # 5      3  (P1 wins on move 3)
        # The formula is (plies + 1) / 2
        turns = (plies + 1) // 2
        
        print(f"Player 1 can force a win in {plies} plies (a ply is one half-move).")
        print(f"The number of turns is calculated from plies as: ({plies} + 1) // 2 = {turns}")
        print(f"The final answer is {turns} turns.")
    elif result == "P2_WIN":
        turns = plies // 2
        print(f"Player 2 can force a win in {turns} turns.")
    else:
        print("The game ends in a draw with optimal play from both sides.")

solve_game()