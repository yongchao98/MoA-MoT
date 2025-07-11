import collections

class GameSolver:
    """
    Solves the described 1D chess-like game to find the shortest forced win for Player 1.
    """

    def __init__(self):
        # Memoization cache to store results for (state, player) tuples.
        self.memo = {}
        # Initial board configuration.
        self.initial_state = {
            'K1': 0, 'N1': 1, 'R1': 2,
            'R2': 5, 'N2': 6, 'K2': 7
        }
        # Player piece sets.
        self.player_pieces = {
            1: ['K1', 'N1', 'R1'],
            2: ['K2', 'N2', 'R2']
        }

    def _is_king_in_check(self, pieces, player):
        """Checks if the specified player's king is under attack by the opponent's rook."""
        king_piece = f'K{player}'
        opponent_rook_piece = f'R{3 - player}'

        if king_piece not in pieces or opponent_rook_piece not in pieces:
            return False

        king_pos = pieces[king_piece]
        rook_pos = pieces[opponent_rook_piece]

        start, end = sorted((king_pos, rook_pos))
        
        # Check for any piece between the king and the rook.
        for piece, pos in pieces.items():
            if piece != king_piece and piece != opponent_rook_piece:
                if start < pos < end:
                    return False  # A piece is blocking the attack.
        
        return True # King is in check.

    def _apply_move(self, pieces, move):
        """Applies a move to a copy of the pieces dictionary and returns the new state."""
        piece_to_move, _, dest_pos = move
        new_pieces = pieces.copy()

        # Check for a capture by finding if any piece is at the destination.
        captured_piece = None
        for p, pos in new_pieces.items():
            if pos == dest_pos:
                captured_piece = p
                break
        
        if captured_piece:
            del new_pieces[captured_piece]
        
        new_pieces[piece_to_move] = dest_pos
        return new_pieces

    def _generate_legal_moves(self, pieces, player):
        """Generates all legal moves for a given player from a given state."""
        legal_moves = []
        occupied_pos = {pos: p for p, pos in pieces.items()}
        
        for piece in self.player_pieces[player]:
            if piece not in pieces:
                continue

            pos = pieces[piece]
            ptype = piece[0]
            
            potential_dests = []
            if ptype == 'K':
                potential_dests.extend([pos - 1, pos + 1])
            elif ptype == 'N':
                potential_dests.extend([pos - 2, pos + 2])
            elif ptype == 'R':
                # Move right
                for d in range(pos + 1, 8):
                    potential_dests.append(d)
                    if d in occupied_pos: break
                # Move left
                for d in range(pos - 1, -1, -1):
                    potential_dests.append(d)
                    if d in occupied_pos: break
            
            for dest in potential_dests:
                # Rule: Check board boundaries
                if not (0 <= dest <= 7):
                    continue
                
                # Rule: Cannot capture own piece
                if dest in occupied_pos and occupied_pos[dest] in self.player_pieces[player]:
                    continue
                
                move = (piece, pos, dest)
                temp_pieces = self._apply_move(pieces, move)

                # Rule: King must not be in check after the move
                if not self._is_king_in_check(temp_pieces, player):
                    legal_moves.append(move)

        return legal_moves

    def _solve_recursive(self, pieces, player):
        """Recursively determines the game outcome from a given state."""
        state_tuple = tuple(sorted(pieces.items()))
        memo_key = (state_tuple, player)

        if memo_key in self.memo:
            return self.memo[memo_key]

        # Base case: A king has been captured
        if f'K{3-player}' not in pieces:
            return ('WIN', 0)
        if f'K{player}' not in pieces:
            return ('LOSS', 0)

        legal_moves = self._generate_legal_moves(pieces, player)

        # Base case: Stalemate (no legal moves) results in a draw
        if not legal_moves:
            return ('DRAW', 0)
        
        # Explore outcomes of all legal moves
        child_outcomes = []
        for move in legal_moves:
            next_pieces = self._apply_move(pieces, move)
            # The outcome for the opponent after this move
            opponent_outcome, opponent_plies = self._solve_recursive(next_pieces, 3 - player)
            child_outcomes.append((opponent_outcome, opponent_plies))
        
        # Determine this state's outcome based on children (minimax logic)
        # Player wants to win. A win is possible if any move leads to an opponent's loss.
        winnable_plies = [1 + plies for outcome, plies in child_outcomes if outcome == 'LOSS']
        if winnable_plies:
            # If win is possible, choose the fastest one.
            result = ('WIN', min(winnable_plies))
        else:
            # If win is not possible, check for a draw.
            is_drawable = any(outcome == 'DRAW' for outcome, plies in child_outcomes)
            if is_drawable:
                result = ('DRAW', float('inf'))
            else:
                # Loss is inevitable. Opponent will try to win, player will try to stall.
                # All moves lead to an opponent win, so this is a loss for the current player.
                # Choose the move that makes the game last as long as possible.
                lossable_plies = [1 + plies for outcome, plies in child_outcomes if outcome == 'WIN']
                result = ('LOSS', max(lossable_plies))

        self.memo[memo_key] = result
        return result

    def solve_game(self):
        """Initiates the solver for the initial game state."""
        # We start as Player 1.
        outcome, plies = self._solve_recursive(self.initial_state, 1)

        if outcome == 'WIN':
            # The problem asks "how many turns", which we interpret as single-player moves (plies).
            print(f"The first player can force a win in {plies} moves.")
        elif outcome == 'DRAW':
            print("The second player can force a draw.")
        else: # LOSS
            print(f"The first player will lose in {plies} moves if the second player plays optimally.")

if __name__ == '__main__':
    solver = GameSolver()
    solver.solve_game()