import sys

def solve_game():
    """
    This function sets up and solves the game, printing the final answer.
    """
    # It's good practice to set a higher recursion limit for deep search problems.
    sys.setrecursionlimit(2000)

    # Memoization table to store results of (board_state, player_turn).
    memo = {}

    # Constants for players and pieces
    P1 = 'P1'
    P2 = 'P2'
    P1_PIECES = {'K1', 'N1', 'R1'}
    P2_PIECES = {'K2', 'N2', 'R2'}

    def get_owner(piece):
        """Determines the owner of a piece."""
        if piece in P1_PIECES:
            return P1
        if piece in P2_PIECES:
            return P2
        return None

    def is_king_in_check(board, player):
        """
        Checks if the specified player's king is under attack by the opponent's rook.
        A king is in check if the opponent's rook has a clear line of sight.
        """
        king_piece = 'K1' if player == P1 else 'K2'
        opponent_rook = 'R2' if player == P1 else 'R1'

        try:
            king_pos = board.index(king_piece)
        except ValueError:
            return False  # King is captured, so not in check.

        try:
            rook_pos = board.index(opponent_rook)
        except ValueError:
            return False  # Opponent's rook is captured, so cannot give check.

        # Check for blocking pieces between the king and the rook.
        start, end = sorted((king_pos, rook_pos))
        for i in range(start + 1, end):
            if board[i] != '':
                return False  # Path is blocked.
        return True

    def generate_legal_moves(board, player):
        """
        Generates all legal moves for a given player.
        A move is legal if it follows piece movement rules and does not leave the
        player's own king in check.
        """
        legal_moves = []
        my_pieces = P1_PIECES if player == P1 else P2_PIECES

        for pos, piece in enumerate(board):
            if piece in my_pieces:
                piece_type = piece[0]
                potential_dests = []
                
                # King: moves one step
                if piece_type == 'K':
                    if pos > 0: potential_dests.append(pos - 1)
                    if pos < 7: potential_dests.append(pos + 1)
                # Knight: moves two steps
                elif piece_type == 'N':
                    if pos > 1: potential_dests.append(pos - 2)
                    if pos < 6: potential_dests.append(pos + 2)
                # Rook: moves any number of unblocked steps
                elif piece_type == 'R':
                    # Move left
                    for dest in range(pos - 1, -1, -1):
                        potential_dests.append(dest)
                        if board[dest] != '': break
                    # Move right
                    for dest in range(pos + 1, 8):
                        potential_dests.append(dest)
                        if board[dest] != '': break

                # Validate each potential destination
                for dest in potential_dests:
                    # A piece cannot move to a square occupied by a friendly piece.
                    if get_owner(board[dest]) == player:
                        continue

                    # Simulate the move to check for king safety.
                    next_board = list(board)
                    next_board[dest] = piece
                    next_board[pos] = ''

                    # A move is illegal if it results in the player's own king being in check.
                    if not is_king_in_check(tuple(next_board), player):
                        legal_moves.append((pos, dest))
        return legal_moves

    def apply_move(board, move):
        """Applies a move to the board and returns the new board state."""
        from_pos, to_pos = move
        piece = board[from_pos]
        next_board = list(board)
        next_board[to_pos] = piece
        next_board[from_pos] = ''
        return tuple(next_board)

    def find_best_outcome(board_tuple, player, path=frozenset()):
        """
        Recursively finds the best possible outcome from a given state using minimax.
        Returns a tuple of (outcome, plies), e.g., ('P1_WIN', 5).
        """
        if (board_tuple, player) in memo:
            return memo[(board_tuple, player)]
        
        if (board_tuple, player) in path:
            return ('DRAW', 0)

        if 'K2' not in board_tuple: return ('P1_WIN', 0)
        if 'K1' not in board_tuple: return ('P2_WIN', 0)

        legal_moves = generate_legal_moves(board_tuple, player)
        if not legal_moves:
            return ('DRAW', 0)

        outcomes = []
        next_player = P2 if player == P1 else P1
        new_path = path | {(board_tuple, player)}
        for move in legal_moves:
            next_board_tuple = apply_move(board_tuple, move)
            outcome, plies = find_best_outcome(next_board_tuple, next_player, new_path)
            outcomes.append((outcome, plies + 1))

        best_result = None
        if player == P1:
            p1_wins = [r for r in outcomes if r[0] == 'P1_WIN']
            if p1_wins:
                best_result = min(p1_wins, key=lambda x: x[1])
            else:
                draws = [r for r in outcomes if r[0] == 'DRAW']
                if draws:
                    best_result = max(draws, key=lambda x: x[1])
                else:
                    best_result = max(outcomes, key=lambda x: x[1])
        else:
            p1_wins = [r for r in outcomes if r[0] == 'P1_WIN']
            if len(p1_wins) == len(outcomes):
                best_result = max(p1_wins, key=lambda x: x[1])
            else:
                non_p1_wins = [r for r in outcomes if r[0] != 'P1_WIN']
                p2_wins = [r for r in non_p1_wins if r[0] == 'P2_WIN']
                if p2_wins:
                    best_result = min(p2_wins, key=lambda x: x[1])
                else:
                    best_result = max(non_p1_wins, key=lambda x: x[1])

        memo[(board_tuple, player)] = best_result
        return best_result

    initial_board = ('K1', 'N1', 'R1', '', '', 'R2', 'N2', 'K2')
    outcome, plies = find_best_outcome(initial_board, P1)

    if outcome == 'P1_WIN':
        # A "turn" for P1 is one of P1's moves.
        # If the game takes M plies and P1 wins, M must be odd.
        # The number of P1 turns is (plies + 1) / 2.
        num_turns = int((plies + 1) / 2)
        print(f"Player 1 can force a win in {num_turns} turns.")
        print(f"<<<{num_turns}>>>")
    else:
        print("Player 1 cannot force a win against an optimal opponent.")

solve_game()