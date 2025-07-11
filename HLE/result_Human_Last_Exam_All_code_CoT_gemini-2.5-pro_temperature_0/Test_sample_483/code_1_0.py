import collections

class GameSolver:
    """
    Solves the 1D chess-like game using a minimax algorithm with memoization.
    """
    def __init__(self):
        # Define pieces and the initial board configuration
        self.P1_PIECES = {"K1", "N1", "R1"}
        self.P2_PIECES = {"K2", "N2", "R2"}
        self.EMPTY = " "
        self.initial_board = ("K1", "N1", "R1", self.EMPTY, self.EMPTY, "R2", "N2", "K2")
        # Memoization table to store results for states already visited
        self.memo = {}

    def is_king_in_check(self, board, player):
        """Checks if the specified player's king is under attack by the opponent's rook."""
        king_piece = "K1" if player == 1 else "K2"
        opp_rook_piece = "R2" if player == 1 else "R1"

        try:
            king_pos = board.index(king_piece)
        except ValueError:
            return False  # King is captured, so not in check

        try:
            rook_pos = board.index(opp_rook_piece)
        except ValueError:
            return False  # Opponent's rook is captured, so no check possible

        # Check if the path between the king and rook is clear
        start, end = sorted((king_pos, rook_pos))
        path = board[start + 1:end]
        return all(p == self.EMPTY for p in path)

    def generate_moves(self, board, player):
        """Generates all legal next board states for the given player."""
        legal_moves = []
        my_pieces = self.P1_PIECES if player == 1 else self.P2_PIECES
        
        for pos, piece in enumerate(board):
            if piece not in my_pieces:
                continue

            piece_type = piece[0]
            
            # Generate potential destination squares based on piece type
            potential_dests = []
            if piece_type == 'K':
                potential_dests.extend([pos - 1, pos + 1])
            elif piece_type == 'N':
                potential_dests.extend([pos - 2, pos + 2])
            elif piece_type == 'R':
                # Scan left for moves
                for dest in range(pos - 1, -1, -1):
                    potential_dests.append(dest)
                    if board[dest] != self.EMPTY:
                        break
                # Scan right for moves
                for dest in range(pos + 1, 8):
                    potential_dests.append(dest)
                    if board[dest] != self.EMPTY:
                        break
            
            # Validate each potential move
            for dest in potential_dests:
                if not (0 <= dest < 8) or board[dest] in my_pieces:
                    continue

                # Create the new board state after the move
                new_board_list = list(board)
                new_board_list[dest] = piece
                new_board_list[pos] = self.EMPTY
                new_board = tuple(new_board_list)

                # A move is legal only if the player's own king is not in check after the move
                if not self.is_king_in_check(new_board, player):
                    legal_moves.append(new_board)
        
        return legal_moves

    def solve(self, board, player):
        """
        Recursively solves the game using minimax.
        Returns the number of turns to a forced win for P1, or float('inf').
        """
        state = (board, player)
        if state in self.memo:
            return self.memo[state]

        # Base cases for recursion termination
        if "K2" not in board:
            return 0  # Player 1 wins, 0 more turns needed
        if "K1" not in board:
            return float('inf') # Player 2 wins, not a P1 winning line

        legal_next_boards = self.generate_moves(board, player)

        if not legal_next_boards:
            # Stalemate: No legal moves, but opponent's king is on board.
            # This is a draw, which is a non-win for P1.
            return float('inf')

        next_player = 2 if player == 1 else 1
        outcomes = [self.solve(b, next_player) for b in legal_next_boards]

        if player == 1:
            # Player 1's turn: minimize the number of turns to win.
            result = 1 + min(outcomes)
        else:  # Player 2's turn
            # Player 2's turn: maximize the number of turns to stall.
            # If P2 can force a draw/win (inf), they will.
            result = 1 + max(outcomes)

        self.memo[state] = result
        return result

    def find_solution(self):
        """Initiates the solver from the starting position."""
        # Start the game with Player 1 to move
        result = self.solve(self.initial_board, 1)
        
        if result == float('inf'):
            print("Player 1 cannot force a win.")
        else:
            print(f"{result}")

if __name__ == "__main__":
    solver = GameSolver()
    solver.find_solution()