import sys

class GameSolver:
    """
    Solves the 1D chess-like game using a minimax algorithm with memoization.
    """
    def __init__(self):
        # Memoization cache for (board_tuple, player) -> (outcome, moves)
        self.memo = {}
        # The game tree can be deep, so a higher recursion limit might be needed.
        # However, the default is usually sufficient for this problem size.
        # sys.setrecursionlimit(3000)

    def get_piece_positions(self, board_tuple):
        """Converts a board tuple to a dictionary of piece positions."""
        positions = {}
        for i, piece in enumerate(board_tuple):
            if piece:
                positions[piece] = i
        return positions

    def is_king_in_check(self, board_tuple, player):
        """Checks if the specified player's king is under attack by the opponent's rook."""
        positions = self.get_piece_positions(board_tuple)
        king_piece = f'K{player}'
        opponent_rook_piece = f'R{2 if player == 1 else 1}'

        if king_piece not in positions or opponent_rook_piece not in positions:
            return False

        king_pos = positions[king_piece]
        rook_pos = positions[opponent_rook_piece]

        start = min(king_pos, rook_pos) + 1
        end = max(king_pos, rook_pos)

        # Check for any blocking pieces between the king and rook
        for i in range(start, end):
            if board_tuple[i] != '':
                return False
        return True

    def generate_legal_moves(self, board_tuple, player):
        """Generates all legal next board states for the given player."""
        board = list(board_tuple)
        legal_boards = set()
        player_suffix = str(player)
        opponent_suffix = str(2 if player == 1 else 1)

        for pos, piece in enumerate(board):
            if not piece or not piece.endswith(player_suffix):
                continue

            piece_type = piece[0]
            
            # Generate moves based on piece type
            potential_moves = []
            if piece_type == 'K': # King
                potential_moves = [pos - 1, pos + 1]
            elif piece_type == 'N': # Knight
                potential_moves = [pos - 2, pos + 2]
            elif piece_type == 'R': # Rook
                # Move right
                for new_pos in range(pos + 1, 8):
                    potential_moves.append(new_pos)
                    if board[new_pos]: break
                # Move left
                for new_pos in range(pos - 1, -1, -1):
                    potential_moves.append(new_pos)
                    if board[new_pos]: break
            
            for new_pos in potential_moves:
                if 0 <= new_pos <= 7:
                    # A move is valid if the destination is empty or an opponent's piece
                    if not board[new_pos] or board[new_pos].endswith(opponent_suffix):
                        next_board_list = list(board)
                        next_board_list[new_pos] = piece
                        next_board_list[pos] = ''
                        next_board_tuple = tuple(next_board_list)
                        
                        # The move is only legal if it doesn't leave the king in check
                        if not self.is_king_in_check(next_board_tuple, player):
                            legal_boards.add(next_board_tuple)
        
        return list(legal_boards)

    def solve(self, board_tuple, player):
        """
        Recursively solves the game using minimax.
        Returns: (outcome, moves_to_outcome)
        Outcome: 1 for P1 win, -1 for P2 win, 0 for draw.
        """
        state_key = (board_tuple, player)
        if state_key in self.memo:
            return self.memo[state_key]

        legal_moves = self.generate_legal_moves(board_tuple, player)

        # Base Case: No legal moves
        if not legal_moves:
            if self.is_king_in_check(board_tuple, player):
                # Player is checkmated, opponent wins in 0 moves from this state.
                result = (-1 if player == 1 else 1, 0)
            else:
                # Player is stalemated. Draw.
                result = (0, 0)
            self.memo[state_key] = result
            return result

        # Recursive Step: Explore outcomes of all legal moves
        opponent = 2 if player == 1 else 1
        outcomes = []
        for next_board in legal_moves:
            outcome, moves = self.solve(next_board, opponent)
            outcomes.append((outcome, moves + 1))

        # Determine the best outcome based on the current player's goal
        if player == 1:  # Player 1 (Maximizer)
            wins = [o for o in outcomes if o[0] == 1]
            if wins: # If there's a winning path, choose the shortest one
                best_result = min(wins, key=lambda x: x[1])
            else:
                draws = [o for o in outcomes if o[0] == 0]
                if draws: # Otherwise, try to force a draw (longest path)
                    best_result = max(draws, key=lambda x: x[1])
                else: # Otherwise, accept loss but prolong it (longest path)
                    best_result = max(outcomes, key=lambda x: x[1])
        else:  # Player 2 (Minimizer)
            wins = [o for o in outcomes if o[0] == -1]
            if wins: # If there's a winning path, choose the shortest one
                best_result = min(wins, key=lambda x: x[1])
            else:
                draws = [o for o in outcomes if o[0] == 0]
                if draws: # Otherwise, try to force a draw (longest path)
                    best_result = max(draws, key=lambda x: x[1])
                else: # Otherwise, P1 wins; prolong the game (longest path)
                    best_result = max(outcomes, key=lambda x: x[1])
        
        self.memo[state_key] = best_result
        return best_result

def main():
    """
    Initializes the game and prints the result.
    """
    # Initial board configuration
    initial_board = ('K1', 'N1', 'R1', '', '', 'R2', 'N2', 'K2')
    
    # Create a solver instance and find the solution
    solver = GameSolver()
    outcome, turns = solver.solve(initial_board, 1)

    # The problem asks for the number of turns for a forced win by Player 1.
    if outcome == 1:
        print(turns)
    else:
        # This case handles scenarios where P1 cannot force a win.
        print("Player 1 cannot force a win.")

if __name__ == '__main__':
    main()