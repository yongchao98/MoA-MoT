import sys

# It's a deep recursion, so we increase the limit.
sys.setrecursionlimit(2000)

class GameSolver:
    """
    Solves the described board game by finding the shortest forced win for Player 1.
    """

    def __init__(self):
        """Initializes the solver with game constants and memoization table."""
        self.memo = {}
        self.INITIAL_BOARD = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
        self.PIECE_OWNERS = {
            'K1': 1, 'N1': 1, 'R1': 1,
            'K2': 2, 'N2': 2, 'R2': 2,
            ' ': None
        }

    def get_opponent(self, player):
        """Returns the opponent of the given player."""
        return 3 - player

    def is_king_in_check(self, board, player):
        """
        Checks if the specified player's King is under attack by the opponent's Rook.
        """
        king_piece = f'K{player}'
        opp_rook_piece = f'R{self.get_opponent(player)}'
        
        try:
            king_pos = board.index(king_piece)
        except ValueError:
            return False  # King is not on the board.

        try:
            opp_rook_pos = board.index(opp_rook_piece)
        except ValueError:
            return False  # Opponent's Rook is not on the board.

        # Check for any blocking pieces between the King and the Rook.
        start = min(king_pos, opp_rook_pos) + 1
        end = max(king_pos, opp_rook_pos)
        
        for i in range(start, end):
            if board[i] != ' ':
                return False  # The path is blocked.
        return True  # The path is clear, so the King is in check.

    def generate_legal_moves(self, board, player):
        """Generates all legal board states reachable in one move by the player."""
        legal_boards = []
        
        for pos, piece in enumerate(board):
            if self.PIECE_OWNERS.get(piece) == player:
                # --- King Moves ---
                if 'K' in piece:
                    for d in [-1, 1]:
                        dest = pos + d
                        if 0 <= dest < 8 and self.PIECE_OWNERS[board[dest]] != player:
                            new_board_list = list(board)
                            new_board_list[dest], new_board_list[pos] = piece, ' '
                            new_board = tuple(new_board_list)
                            if not self.is_king_in_check(new_board, player):
                                legal_boards.append(new_board)
                # --- Knight Moves ---
                elif 'N' in piece:
                    for d in [-2, 2]:
                        dest = pos + d
                        if 0 <= dest < 8 and self.PIECE_OWNERS[board[dest]] != player:
                            new_board_list = list(board)
                            new_board_list[dest], new_board_list[pos] = piece, ' '
                            new_board = tuple(new_board_list)
                            if not self.is_king_in_check(new_board, player):
                                legal_boards.append(new_board)
                # --- Rook Moves ---
                elif 'R' in piece:
                    # Move Left
                    for dest in range(pos - 1, -1, -1):
                        target_owner = self.PIECE_OWNERS[board[dest]]
                        if target_owner == player: break
                        new_board_list = list(board)
                        new_board_list[dest], new_board_list[pos] = piece, ' '
                        new_board = tuple(new_board_list)
                        if not self.is_king_in_check(new_board, player):
                            legal_boards.append(new_board)
                        if target_owner is not None: break
                    # Move Right
                    for dest in range(pos + 1, 8):
                        target_owner = self.PIECE_OWNERS[board[dest]]
                        if target_owner == player: break
                        new_board_list = list(board)
                        new_board_list[dest], new_board_list[pos] = piece, ' '
                        new_board = tuple(new_board_list)
                        if not self.is_king_in_check(new_board, player):
                            legal_boards.append(new_board)
                        if target_owner is not None: break
        return legal_boards

    def minimax(self, state):
        """
        Recursively determines the game outcome from the given state using minimax.
        Returns: (winner, turns)
        """
        board, player = state
        state_key = (board, player)

        if state_key in self.memo:
            return self.memo[state_key]

        # --- Base Cases ---
        opponent = self.get_opponent(player)
        if f'K{opponent}' not in board: return (player, 0)
        if f'K{player}' not in board: return (opponent, 0)
        
        legal_next_boards = self.generate_legal_moves(board, player)

        if not legal_next_boards:
            return (0, 0) # Stalemate results in a draw

        # --- Recursive Step ---
        results = []
        for next_board in legal_next_boards:
            outcome, turns = self.minimax((next_board, opponent))
            results.append((outcome, turns + 1))

        # Determine the best outcome based on the current player's goals
        wins = [r for r in results if r[0] == player]
        draws = [r for r in results if r[0] == 0]
        losses = [r for r in results if r[0] == opponent]

        # Player wants to win fast
        if wins:
            best_res = min(wins, key=lambda x: x[1])
        # If winning is not possible, player settles for a draw
        elif draws:
            best_res = (0, 0)
        # If losing is inevitable, player wants to stall as long as possible
        else:
            best_res = max(losses, key=lambda x: x[1])
            
        self.memo[state_key] = best_res
        return best_res

    def solve_game(self):
        """Solves the game from the initial state and prints the result."""
        initial_state = (self.INITIAL_BOARD, 1)
        outcome, turns = self.minimax(initial_state)

        if outcome == 1:
            print(f"Result: Player 1 can force a win.")
            print(f"The number of turns for Player 1 to force a win is: {turns}")
            # This is the final answer for the user's question
            self.final_answer = turns
        elif outcome == 2:
            print(f"Result: Player 2 can force a win in {turns} turns.")
            self.final_answer = -1 # Indicates P1 cannot force a win
        else:
            print("Result: The game is a forced draw.")
            self.final_answer = 0 # Indicates a draw


if __name__ == '__main__':
    solver = GameSolver()
    solver.solve_game()
