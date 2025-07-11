import functools

class GameSolver:
    """
    Solves the described board game to find the shortest forced win for Player 1.
    """
    PIECES = {
        'K1': 1, 'N1': 1, 'R1': 1,
        'K2': 2, 'N2': 2, 'R2': 2
    }

    def __init__(self):
        # Apply memoization to the main solver function to cache results of game states.
        # This is crucial for performance.
        self.solve = functools.lru_cache(maxsize=None)(self._solve)

    def get_piece_positions(self, board):
        """Helper function to get a map of pieces to their positions."""
        positions = {}
        for i, piece in enumerate(board):
            if piece:
                positions[piece] = i
        return positions

    def is_king_in_check(self, board, player):
        """
        Checks if the specified player's King is under attack by the opponent's Rook.
        A King is in check if the Rook has a clear line of sight.
        """
        positions = self.get_piece_positions(board)
        king_piece = f'K{player}'
        opponent_rook_piece = f'R{3 - player}'

        if king_piece not in positions or opponent_rook_piece not in positions:
            return False  # King or attacking rook is off the board.

        king_pos = positions[king_piece]
        rook_pos = positions[opponent_rook_piece]

        start, end = min(king_pos, rook_pos), max(king_pos, rook_pos)
        # Check for any blocking pieces between the King and the Rook.
        for i in range(start + 1, end):
            if board[i] is not None:
                return False  # Path is blocked.
        return True # Path is clear, King is in check.

    def generate_moves(self, board, player):
        """Generates all possible moves for a player, without checking for king safety yet."""
        moves = []
        player_pieces_on_board = [(p, i) for i, p in enumerate(board) if p and self.PIECES.get(p) == player]

        for piece, pos in player_pieces_on_board:
            # King: moves 1 step
            if 'K' in piece:
                for d in [-1, 1]:
                    new_pos = pos + d
                    if 0 <= new_pos <= 7 and (board[new_pos] is None or self.PIECES.get(board[new_pos]) != player):
                        new_board = list(board)
                        new_board[pos], new_board[new_pos] = None, piece
                        moves.append(tuple(new_board))
            # Knight: moves 2 steps
            elif 'N' in piece:
                for d in [-2, 2]:
                    new_pos = pos + d
                    if 0 <= new_pos <= 7 and (board[new_pos] is None or self.PIECES.get(board[new_pos]) != player):
                        new_board = list(board)
                        new_board[pos], new_board[new_pos] = None, piece
                        moves.append(tuple(new_board))
            # Rook: moves any number of steps until blocked
            elif 'R' in piece:
                for direction in [-1, 1]:
                    for i in range(1, 8):
                        new_pos = pos + direction * i
                        if not (0 <= new_pos <= 7): break
                        
                        target = board[new_pos]
                        if target and self.PIECES.get(target) == player: break
                        
                        new_board = list(board)
                        new_board[pos], new_board[new_pos] = None, piece
                        moves.append(tuple(new_board))
                        
                        if target and self.PIECES.get(target) != player: break
        return moves

    def get_legal_moves(self, board, player):
        """Filters generated moves to only include legal ones (King is not left in check)."""
        return [move for move in self.generate_moves(board, player) if not self.is_king_in_check(move, player)]

    def _solve(self, board, player):
        """
        The core minimax solver.
        Returns a tuple (outcome, plies):
        - outcome: 1 for P1 win, 2 for P2 win, 0 for draw.
        - plies: The number of moves (plies) to reach the outcome.
        """
        # Base case: A king has been captured.
        positions = self.get_piece_positions(board)
        if 'K2' not in positions: return (1, 0)
        if 'K1' not in positions: return (2, 0)

        # Base case: No legal moves available.
        legal_moves = self.get_legal_moves(board, player)
        if not legal_moves:
            # Checkmate if in check, otherwise stalemate.
            if self.is_king_in_check(board, player):
                return (3 - player, 0)  # Current player is checkmated and loses.
            else:
                return (0, float('inf'))  # Stalemate (draw).

        # Recursive step: evaluate all legal moves.
        outcomes = [self.solve(move, 3 - player) for move in legal_moves]

        # Player 1's turn: try to win as fast as possible.
        if player == 1:
            wins = [(o, p) for o, p in outcomes if o == 1]
            if wins:
                min_plies = min(p for o, p in wins)
                return (1, 1 + min_plies) # Choose the fastest win.
            
            if any(o == 0 for o, p in outcomes):
                return (0, float('inf')) # Can force a draw.
            
            # Must lose, so choose the move that prolongs the game (P2 wins slowest).
            max_plies = max(p for o, p in outcomes)
            return (2, 1 + max_plies)

        # Player 2's turn: try to win, or stall a loss as long as possible.
        else:
            wins = [(o, p) for o, p in outcomes if o == 2]
            if wins:
                min_plies = min(p for o, p in wins)
                return (2, 1 + min_plies) # Choose the fastest win.

            if any(o == 0 for o, p in outcomes):
                return (0, float('inf')) # Can force a draw.
            
            # Must lose, so choose the move that stalls the longest (P1 wins slowest).
            max_plies = max(p for o, p in outcomes)
            return (1, 1 + max_plies)

    def find_solution(self):
        """Initializes the game and prints the final result."""
        initial_board = ('K1', 'N1', 'R1', None, None, 'R2', 'N2', 'K2')
        initial_player = 1
        
        outcome, plies = self.solve(initial_board, initial_player)

        if outcome == 1:
            print(plies)
        elif outcome == 2:
            print(f"Player 2 can force a win in {plies} moves.")
        else:
            print("The game results in a forced draw.")

if __name__ == '__main__':
    solver = GameSolver()
    solver.find_solution()