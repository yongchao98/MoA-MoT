import functools

class GameSolver:
    """
    A class to solve the described board game.
    """
    # Constants for players and game pieces
    P1 = 1
    P2 = 2
    DRAW = 0
    INITIAL_BOARD = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
    PIECE_OWNERS = {
        'K1': P1, 'N1': P1, 'R1': P1,
        'K2': P2, 'N2': P2, 'R2': P2
    }

    def get_opponent(self, player):
        """Returns the opponent of the given player."""
        return self.P2 if player == self.P1 else self.P1

    def is_in_check(self, board, player):
        """Checks if the given player's King is under attack by the opponent's Rook."""
        king_piece = 'K1' if player == self.P1 else 'K2'
        opponent_rook_piece = 'R2' if player == self.P1 else 'R1'

        try:
            king_pos = board.index(king_piece)
        except ValueError:
            return False  # King is not on the board, so not in check.

        try:
            rook_pos = board.index(opponent_rook_piece)
        except ValueError:
            return False  # Opponent's Rook is not on the board.

        # Check for a clear path between the King and the Rook.
        start, end = sorted((king_pos, rook_pos))
        for i in range(start + 1, end):
            if board[i] != ' ':
                return False  # Path is blocked.
        
        return True # Path is clear, King is in check.

    def generate_legal_moves(self, board, player):
        """Generates all legal next board states for the given player."""
        legal_boards = []
        for pos, piece in enumerate(board):
            if piece == ' ' or self.PIECE_OWNERS.get(piece) != player:
                continue

            piece_type = piece[0]
            potential_dests = []

            if piece_type == 'K':
                potential_dests = [pos - 1, pos + 1]
            elif piece_type == 'N':
                potential_dests = [pos - 2, pos + 2]
            elif piece_type == 'R':
                for i in range(pos - 1, -1, -1):
                    potential_dests.append(i)
                    if board[i] != ' ': break
                for i in range(pos + 1, 8):
                    potential_dests.append(i)
                    if board[i] != ' ': break
            
            for dest in potential_dests:
                if not (0 <= dest < 8): continue
                if board[dest] != ' ' and self.PIECE_OWNERS.get(board[dest]) == player: continue
                
                next_board_list = list(board)
                next_board_list[dest], next_board_list[pos] = piece, ' '
                next_board = tuple(next_board_list)
                
                if not self.is_in_check(next_board, player):
                    legal_boards.append(next_board)
                    
        return legal_boards

    @functools.lru_cache(maxsize=None)
    def solve(self, board, player, history):
        """
        Recursively solves the game state using minimax, returning (winner, ply_to_end).
        - winner: P1, P2, or DRAW
        - ply_to_end: Number of moves (ply) until the game ends.
        - history: A frozenset of visited states to detect draws by repetition.
        """
        state = (board, player)
        if state in history:
            return (self.DRAW, 0)

        legal_next_boards = self.generate_legal_moves(board, player)

        if not legal_next_boards:
            if self.is_in_check(board, player):
                return (self.get_opponent(player), 0) # Checkmate
            else:
                return (self.DRAW, 0) # Stalemate

        new_history = history | {state}
        outcomes = [self.solve(b, self.get_opponent(player), new_history) for b in legal_next_boards]

        if player == self.P1: # Player 1 (maximizer) wants to win
            p1_wins = [(w, p) for w, p in outcomes if w == self.P1]
            if p1_wins: # Find the fastest win
                min_ply = min(p for _, p in p1_wins)
                result = (self.P1, min_ply + 1)
            else: # If no win is possible, try to draw
                draws = [(w, p) for w, p in outcomes if w == self.DRAW]
                if draws:
                    max_ply = max(p for _, p in draws)
                    result = (self.DRAW, max_ply + 1)
                else: # Must lose, so stall as long as possible
                    max_ply = max(p for _, p in outcomes)
                    result = (self.P2, max_ply + 1)
        else: # Player 2 (minimizer) wants to win/draw/stall
            p2_wins = [(w, p) for w, p in outcomes if w == self.P2]
            if p2_wins: # Find the fastest win
                min_ply = min(p for _, p in p2_wins)
                result = (self.P2, min_ply + 1)
            else:
                draws = [(w, p) for w, p in outcomes if w == self.DRAW]
                if draws:
                    max_ply = max(p for _, p in draws)
                    result = (self.DRAW, max_ply + 1)
                else: # Must lose to P1, so stall as long as possible
                    max_ply = max(p for _, p in outcomes)
                    result = (self.P1, max_ply + 1)

        return result

    def find_solution(self):
        """
        Starts the solver from the initial game state.
        """
        winner, ply = self.solve(self.INITIAL_BOARD, self.P1, frozenset())
        
        if winner == self.P1:
            # A "turn" for Player 1 consists of one of their moves.
            # ply=1 -> P1 moves once -> 1 turn.
            # ply=3 -> P1,P2,P1 -> 2 turns for P1.
            # Formula for P1 turns is (ply + 1) // 2
            p1_turns = (ply + 1) // 2
            print(p1_turns)
        elif winner == self.P2:
            print("Player 2 can force a win.")
        else:
            print("The game results in a forced draw.")

if __name__ == '__main__':
    solver = GameSolver()
    solver.find_solution()