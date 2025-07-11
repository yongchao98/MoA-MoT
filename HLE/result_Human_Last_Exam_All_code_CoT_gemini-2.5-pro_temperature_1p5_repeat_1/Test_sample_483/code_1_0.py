import sys

# It's possible the recursion depth will be high for some game lines,
# so we increase the limit.
sys.setrecursionlimit(2000)

class LinearChessSolver:
    """
    Solves the described linear chess puzzle by finding the shortest forced win for Player 1
    using a minimax search algorithm with memoization.
    """

    def __init__(self):
        # Memoization cache to store results of (board, player) states
        self.memo = {}
        # Constants for players and pieces
        self.PLAYER_1 = 'P1'
        self.PLAYER_2 = 'P2'
        self.PIECES = {
            'P1': {'K1', 'N1', 'R1'},
            'P2': {'K2', 'N2', 'R2'}
        }
        self.KINGS = {'P1': 'K1', 'P2': 'K2'}
        self.ROOKS = {'P1': 'R1', 'P2': 'R2'}

    def get_opponent(self, player):
        """Returns the opponent of the given player."""
        return self.PLAYER_2 if player == self.PLAYER_1 else self.PLAYER_1

    def is_in_check(self, board, player):
        """
        Checks if the specified player's King is under attack by the opponent's Rook.
        """
        king_piece = self.KINGS[player]
        try:
            king_pos = board.index(king_piece)
        except ValueError:
            return False  # King is captured, so not considered "in check".

        opponent = self.get_opponent(player)
        opponent_rook_piece = self.ROOKS[opponent]
        try:
            opponent_rook_pos = board.index(opponent_rook_piece)
        except ValueError:
            return False  # Opponent's Rook is captured.

        # Check for a clear line of sight between the King and the Rook
        start = min(king_pos, opponent_rook_pos)
        end = max(king_pos, opponent_rook_pos)
        for i in range(start + 1, end):
            if board[i] != ' ':
                return False  # Path is blocked by another piece.
        return True

    def generate_legal_moves(self, board, player):
        """
        Generates all legal moves for the given player.
        A move is legal if it follows piece rules and does not result in the player's
        own King being in check.
        """
        legal_moves = []
        player_pieces = self.PIECES[player]

        for pos, piece in enumerate(board):
            if piece not in player_pieces:
                continue

            targets = []
            if 'K' in piece: # King moves one step
                targets.extend([pos - 1, pos + 1])
            elif 'N' in piece: # Knight moves two steps
                targets.extend([pos - 2, pos + 2])
            elif 'R' in piece: # Rook moves any number of steps in a clear line
                # Scan right
                for i in range(pos + 1, 8):
                    targets.append(i)
                    if board[i] != ' ': break
                # Scan left
                for i in range(pos - 1, -1, -1):
                    targets.append(i)
                    if board[i] != ' ': break
            
            # Check validity and king safety of each potential move
            for new_pos in targets:
                if 0 <= new_pos < 8 and board[new_pos] not in player_pieces:
                    # Create a temporary board to check for king safety
                    temp_board = list(board)
                    temp_board[new_pos] = piece
                    temp_board[pos] = ' '
                    if not self.is_in_check(tuple(temp_board), player):
                        legal_moves.append((pos, new_pos))
        return legal_moves

    def apply_move(self, board, move):
        """Applies a move to the board and returns the new board state."""
        from_pos, to_pos = move
        new_board_list = list(board)
        piece = new_board_list[from_pos]
        new_board_list[to_pos] = piece
        new_board_list[from_pos] = ' '
        return tuple(new_board_list)

    def solve(self, board, player):
        """
        Recursively solves the game from a given state using minimax.
        Returns a tuple of (outcome, plies).
        """
        state = (board, player)
        if state in self.memo:
            return self.memo[state]

        # Base Cases: Win/Loss by capture or Stalemate
        if self.KINGS[self.PLAYER_2] not in board: return ('P1_WIN', 0)
        if self.KINGS[self.PLAYER_1] not in board: return ('P2_WIN', 0)
        
        moves = self.generate_legal_moves(board, player)
        if not moves: return ('DRAW', 0)

        # Recursive Step: Explore all moves
        opponent = self.get_opponent(player)
        outcomes = []
        for move in moves:
            new_board = self.apply_move(board, move)
            result, plies = self.solve(new_board, opponent)
            outcomes.append((result, plies + 1))

        # Minimax logic to determine the best outcome
        if player == self.PLAYER_1:  # Player 1 (maximizer)
            p1_wins = [p for r, p in outcomes if r == 'P1_WIN']
            if p1_wins: return ('P1_WIN', min(p1_wins)) # Choose the fastest win

            draws = [p for r, p in outcomes if r == 'DRAW']
            if draws: return ('DRAW', max(draws)) # P2 will prolong a draw

            p2_wins = [p for r, p in outcomes if r == 'P2_WIN']
            return ('P2_WIN', max(p2_wins)) # P1 must lose, so stalls as long as possible
        else:  # Player 2 (minimizer)
            p2_wins = [p for r, p in outcomes if r == 'P2_WIN']
            if p2_wins: return ('P2_WIN', min(p2_wins)) # Choose the fastest win

            draws = [p for r, p in outcomes if r == 'DRAW']
            if draws: return ('DRAW', max(draws)) # P2 can force a draw, will prolong it
            
            p1_wins = [p for r, p in outcomes if r == 'P1_WIN']
            return ('P1_WIN', max(p1_wins)) # P2 must lose, so stalls as long as possible

    def find_solution(self):
        """
        Sets up the initial board, starts the solver, and prints the final answer.
        """
        initial_board = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
        outcome, plies = self.solve(initial_board, self.PLAYER_1)
        
        if outcome == 'P1_WIN':
            # A "turn" for Player 1 to win consists of N pairs of moves plus one final move.
            # 1 ply (P1 wins on move 1) -> 1st turn.
            # 3 plies (P1, P2, P1 wins) -> 2nd turn.
            # Turns = ceil(plies / 2)
            num_turns = (plies + 1) // 2
            print(num_turns)
        else:
            print(f"Player 1 cannot force a win. Optimal play leads to a {outcome}.")

if __name__ == '__main__':
    solver = LinearChessSolver()
    solver.find_solution()