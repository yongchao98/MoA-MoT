import copy

class Connect4Solver:
    """
    A class to solve a given Connect 4 board state.
    """
    def __init__(self, board_state, player='O'):
        self.rows = 6
        self.cols = 7
        self.board = board_state
        self.player = player
        self.opponent = 'X' if player == 'O' else 'O'

    def get_next_open_row(self, board, col):
        """Finds the next empty row in a given column."""
        for r in range(self.rows - 1, -1, -1):
            if board[r][col] == '.':
                return r
        return None

    def get_valid_moves(self, board):
        """Gets a list of columns that are not full."""
        return [c for c in range(self.cols) if self.get_next_open_row(board, c) is not None]

    def check_win(self, board, player):
        """Checks if the specified player has won the game."""
        # Check horizontal locations for win
        for c in range(self.cols - 3):
            for r in range(self.rows):
                if all(board[r][c + i] == player for i in range(4)):
                    return True

        # Check vertical locations for win
        for c in range(self.cols):
            for r in range(self.rows - 3):
                if all(board[r + i][c] == player for i in range(4)):
                    return True

        # Check positively sloped diagonals
        for c in range(self.cols - 3):
            for r in range(self.rows - 3):
                if all(board[r + i][c + i] == player for i in range(4)):
                    return True

        # Check negatively sloped diagonals
        for c in range(self.cols - 3):
            for r in range(3, self.rows):
                if all(board[r - i][c + i] == player for i in range(4)):
                    return True
        return False

    def find_fastest_win(self):
        """
        Finds moves that guarantee a win on the next turn.
        A move is optimal if it creates two or more immediate winning threats
        that the opponent cannot block in a single move.
        """
        optimal_moves = []
        
        # 1. Get all possible moves for the current player
        possible_moves = self.get_valid_moves(self.board)
        
        for move_col in possible_moves:
            # 2. Simulate making the move
            temp_board = copy.deepcopy(self.board)
            move_row = self.get_next_open_row(temp_board, move_col)
            temp_board[move_row][move_col] = self.player

            # If this move itself wins, it's the best. But the problem implies a win in 1 turn.
            if self.check_win(temp_board, self.player):
                # This is a win-in-0, which is the fastest.
                # In this problem, no such move exists, but we handle it for completeness.
                return [(move_col, move_row)]

            # 3. Count winning follow-up moves (threats)
            winning_followups = set()
            possible_next_moves = self.get_valid_moves(temp_board)
            
            for followup_col in possible_next_moves:
                final_board = copy.deepcopy(temp_board)
                final_row = self.get_next_open_row(final_board, followup_col)
                final_board[final_row][followup_col] = self.player
                
                if self.check_win(final_board, self.player):
                    winning_followups.add(followup_col)

            # 4. If there are >= 2 winning follow-ups, it's a guaranteed win
            # We must also check that our move doesn't allow the opponent to win immediately.
            if len(winning_followups) >= 2:
                opponent_can_win = False
                opponent_moves = self.get_valid_moves(temp_board)
                for opp_move_col in opponent_moves:
                    opp_board = copy.deepcopy(temp_board)
                    opp_move_row = self.get_next_open_row(opp_board, opp_move_col)
                    opp_board[opp_move_row][opp_move_col] = self.opponent
                    if self.check_win(opp_board, self.opponent):
                        opponent_can_win = True
                        break
                
                if not opponent_can_win:
                    optimal_moves.append((move_col, move_row))
        
        return optimal_moves


# --- Main Execution ---

# Board representation from the problem (top-down, 0-indexed)
# Row 0 is row '1' in the diagram, Row 5 is row '6'.
# Col 0 is 'a', Col 6 is 'g'.
initial_board = [
    ['.', '.', '.', '.', '.', '.', '.'],  # Row 1
    ['.', '.', '.', '.', '.', '.', '.'],  # Row 2
    ['.', '.', '.', '.', '.', '.', '.'],  # Row 3
    ['.', '.', '.', 'O', 'O', '.', '.'],  # Row 4
    ['O', '.', 'X', 'O', 'X', 'X', 'X'],  # Row 5
    ['X', 'O', 'O', 'X', 'X', 'O', 'X']   # Row 6
]

solver = Connect4Solver(initial_board, 'O')
winning_moves = solver.find_fastest_win()

# Format the output as per the user's request
if not winning_moves:
    print("No optimal move found to win in the next turn.")
else:
    col_map = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g'}
    result_coords = []
    for col_idx, row_idx in winning_moves:
        # Convert board row index (0-5) to game row number (1-6)
        row_num = 6 - row_idx
        result_coords.append(f"{col_map[col_idx]}{row_num}")
    
    print(", ".join(result_coords))