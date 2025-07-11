import sys

# Disable print to stderr for the final output
# sys.stderr = open(os.devnull, 'w')

class Connect4Solver:
    """
    A class to solve a specific Connect 4 board state.
    """
    def __init__(self, board_state, player):
        """
        Initializes the solver with a board state.
        Args:
            board_state (list of list of str): The 6x7 game board.
            player (str): The player whose turn it is ('O' or 'X').
        """
        self.rows = 6
        self.cols = 7
        self.board = board_state
        self.player = player
        self.opponent = 'X' if player == 'O' else 'O'

    def get_valid_moves(self, board):
        """
        Gets a list of column indices where a piece can be legally dropped.
        """
        return [c for c in range(self.cols) if board[0][c] == '.']

    def get_landing_row(self, board, col):
        """
        Finds the row a piece will land in for a given column.
        Returns -1 if the column is full.
        """
        for r in range(self.rows - 1, -1, -1):
            if board[r][col] == '.':
                return r
        return -1

    def make_move(self, board, col, player):
        """
        Creates a new board state after a player makes a move.
        Returns None if the move is invalid.
        """
        row = self.get_landing_row(board, col)
        if row == -1:
            return None
        new_board = [r[:] for r in board]
        new_board[row][col] = player
        return new_board

    def check_win(self, board, player):
        """
        Checks if the specified player has won on the given board.
        """
        # Horizontal check
        for r in range(self.rows):
            for c in range(self.cols - 3):
                if all(board[r][c + i] == player for i in range(4)):
                    return True
        # Vertical check
        for c in range(self.cols):
            for r in range(self.rows - 3):
                if all(board[r + i][c] == player for i in range(4)):
                    return True
        # Diagonal (top-left to bottom-right)
        for r in range(self.rows - 3):
            for c in range(self.cols - 3):
                if all(board[r + i][c + i] == player for i in range(4)):
                    return True
        # Diagonal (bottom-left to top-right)
        for r in range(3, self.rows):
            for c in range(self.cols - 3):
                if all(board[r - i][c + i] == player for i in range(4)):
                    return True
        return False

    def find_winning_moves(self, board, player):
        """
        Finds all column moves that result in an immediate win for the player.
        """
        winning_moves = []
        for col in self.get_valid_moves(board):
            temp_board = self.make_move(board, col, player)
            if self.check_win(temp_board, player):
                winning_moves.append(col)
        return winning_moves

    def solve(self):
        """
        Finds all optimal moves for the current player to win as fast as possible.
        """
        optimal_moves_cols = set()
        possible_moves = self.get_valid_moves(self.board)

        for move_col in possible_moves:
            # Simulate the player's move
            board_after_move = self.make_move(self.board, move_col, self.player)

            # A move is bad if it allows the opponent to win immediately.
            if self.find_winning_moves(board_after_move, self.opponent):
                continue

            # Find all winning responses for the player on their *next* turn.
            player_winning_responses = self.find_winning_moves(board_after_move, self.player)

            # Case 1: Simple Fork. The move creates two or more winning threats.
            # The opponent can only block one.
            if len(player_winning_responses) >= 2:
                optimal_moves_cols.add(move_col)
                continue

            # Case 2: Forced Move Trap. The move creates one threat, forcing a block
            # that enables a different winning move for the player.
            if len(player_winning_responses) == 1:
                forced_block_col = player_winning_responses[0]
                
                # Assume opponent must block. They have no winning move of their own.
                # Simulate the opponent's forced block.
                board_after_block = self.make_move(board_after_move, forced_block_col, self.opponent)

                # Check if the player has a new winning move after the block.
                if self.find_winning_moves(board_after_block, self.player):
                    optimal_moves_cols.add(move_col)

        # Convert column indices to the required "a1" format.
        col_map = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        optimal_moves_str = []
        for col in sorted(list(optimal_moves_cols)):
            row_idx = self.get_landing_row(self.board, col)
            row_num = self.rows - row_idx
            optimal_moves_str.append(f"{col_map[col]}{row_num}")
        
        print(", ".join(optimal_moves_str))


# The board state from the problem description
# '.' represents an empty cell.
initial_board = [
    ['.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', 'O', '.', '.'],
    ['O', '.', 'X', 'O', 'X', 'X', 'X'],
    ['X', 'O', 'O', 'X', 'X', 'O', 'X']
]

# Create solver and find the optimal moves for player 'O'
solver = Connect4Solver(initial_board, 'O')
solver.solve()