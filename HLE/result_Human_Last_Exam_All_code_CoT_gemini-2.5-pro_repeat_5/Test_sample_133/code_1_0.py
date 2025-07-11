import collections

class Connect4Solver:
    def __init__(self, initial_board):
        self.initial_board = initial_board
        self.rows = 6
        self.cols = 7
        self.player = 'O'
        self.opponent = 'X'

    def get_valid_moves(self, board):
        """Returns a list of column indices that are not full."""
        return [c for c in range(self.cols) if board[0][c] == '.']

    def get_next_open_row(self, board, col):
        """Returns the row index where a piece would land in the given column."""
        for r in range(self.rows - 1, -1, -1):
            if board[r][col] == '.':
                return r
        return -1

    def make_move(self, board, col, player):
        """Returns a new board state after a move is made."""
        if col is None or self.get_next_open_row(board, col) == -1:
            return None
        temp_board = [row[:] for row in board]
        row = self.get_next_open_row(temp_board, col)
        temp_board[row][col] = player
        return temp_board

    def check_win(self, board, player):
        """Checks if the given player has won on the board."""
        # Check horizontal
        for r in range(self.rows):
            for c in range(self.cols - 3):
                if all(board[r][c + i] == player for i in range(4)):
                    return True
        # Check vertical
        for r in range(self.rows - 3):
            for c in range(self.cols):
                if all(board[r + i][c] == player for i in range(4)):
                    return True
        # Check diagonal (down-right)
        for r in range(self.rows - 3):
            for c in range(self.cols - 3):
                if all(board[r + i][c + i] == player for i in range(4)):
                    return True
        # Check diagonal (up-right)
        for r in range(3, self.rows):
            for c in range(self.cols - 3):
                if all(board[r - i][c + i] == player for i in range(4)):
                    return True
        return False
    
    def find_winning_spots(self, board, player):
        """Finds all empty spots that would result in a win for the player."""
        winning_spots = []
        for c in range(self.cols):
            r = self.get_next_open_row(board, c)
            if r != -1:
                temp_board = [row[:] for row in board]
                temp_board[r][c] = player
                if self.check_win(temp_board, player):
                    winning_spots.append((r,c))
        return winning_spots

    def find_optimal_moves(self):
        """Finds all moves that lead to a win in the shortest number of turns."""
        valid_moves = self.get_valid_moves(self.initial_board)
        
        # Level 1: Check for immediate winning moves
        immediate_wins = []
        for move_col in valid_moves:
            board_after_move = self.make_move(self.initial_board, move_col, self.player)
            if self.check_win(board_after_move, self.player):
                immediate_wins.append(move_col)
        
        if immediate_wins:
            return self.format_moves(immediate_wins)

        # Level 2: Check for moves that guarantee a win on the next turn
        guaranteed_wins = []
        for move_col in valid_moves:
            board_after_O_move = self.make_move(self.initial_board, move_col, self.player)
            
            # This move is a guaranteed win if for all opponent responses, we still have a winning move.
            is_guaranteed = True
            opponent_moves = self.get_valid_moves(board_after_O_move)
            if not opponent_moves: # Opponent has no moves left, it's a draw.
                is_guaranteed = False

            for opp_move_col in opponent_moves:
                board_after_X_move = self.make_move(board_after_O_move, opp_move_col, self.opponent)
                if not self.find_winning_spots(board_after_X_move, self.player):
                    is_guaranteed = False
                    break # Opponent found a block
            
            if is_guaranteed:
                guaranteed_wins.append(move_col)

        return self.format_moves(guaranteed_wins)

    def format_moves(self, cols):
        """Converts column indices to algebraic notation like 'c4'."""
        move_names = []
        for c in cols:
            r_idx = self.get_next_open_row(self.initial_board, c)
            row_num = self.rows - r_idx
            col_char = chr(ord('a') + c)
            move_names.append(f"{col_char}{row_num}")
        return sorted(move_names)

# --- Main execution ---
# Board representation (row 0 is the top row)
board_state = [
  ['.', '.', '.', '.', '.', '.', '.'],
  ['.', '.', '.', '.', '.', '.', '.'],
  ['.', '.', '.', '.', '.', '.', '.'],
  ['.', '.', '.', 'O', 'O', '.', '.'],
  ['O', '.', 'X', 'O', 'X', 'X', 'X'],
  ['X', 'O', 'O', 'X', 'X', 'O', 'X']
]

solver = Connect4Solver(board_state)
optimal_moves = solver.find_optimal_moves()

print(", ".join(optimal_moves))