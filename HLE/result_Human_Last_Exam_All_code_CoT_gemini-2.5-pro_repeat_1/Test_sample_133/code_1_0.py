import copy

class Connect4Solver:
    def __init__(self, board_state):
        self.board = board_state
        self.ROWS = 6
        self.COLS = 7
        self.PLAYER_X = 1
        self.PLAYER_O = 2

    def get_next_open_row(self, board, col):
        """Gets the next open row for a given column."""
        for r in range(self.ROWS - 1, -1, -1):
            if board[r][col] == 0:
                return r
        return -1

    def is_valid_location(self, board, col):
        """Checks if a column is not full."""
        return board[0][col] == 0

    def get_valid_locations(self, board):
        """Gets a list of all valid columns to play in."""
        valid_locations = []
        for col in range(self.COLS):
            if self.is_valid_location(board, col):
                valid_locations.append(col)
        return valid_locations

    def check_win(self, board, player):
        """Checks if the given player has won."""
        # Check horizontal locations for win
        for c in range(self.COLS - 3):
            for r in range(self.ROWS):
                if board[r][c] == player and board[r][c+1] == player and board[r][c+2] == player and board[r][c+3] == player:
                    return True

        # Check vertical locations for win
        for c in range(self.COLS):
            for r in range(self.ROWS - 3):
                if board[r][c] == player and board[r+1][c] == player and board[r+2][c] == player and board[r+3][c] == player:
                    return True

        # Check positively sloped diaganols
        for c in range(self.COLS - 3):
            for r in range(self.ROWS - 3):
                if board[r][c] == player and board[r+1][c+1] == player and board[r+2][c+2] == player and board[r+3][c+3] == player:
                    return True

        # Check negatively sloped diaganols
        for c in range(self.COLS - 3):
            for r in range(3, self.ROWS):
                if board[r][c] == player and board[r-1][c+1] == player and board[r-2][c+2] == player and board[r-3][c+3] == player:
                    return True
        return False

    def solve(self):
        """Finds the optimal moves for player O."""
        # --- Check for immediate wins (Mate in 1) ---
        winning_moves_m1 = []
        possible_moves = self.get_valid_locations(self.board)
        for col in possible_moves:
            temp_board = copy.deepcopy(self.board)
            row = self.get_next_open_row(temp_board, col)
            temp_board[row][col] = self.PLAYER_O
            if self.check_win(temp_board, self.PLAYER_O):
                move_str = f"{chr(ord('a') + col)}{self.ROWS - row}"
                winning_moves_m1.append(move_str)
        
        if winning_moves_m1:
            print(", ".join(sorted(winning_moves_m1)))
            return

        # --- Check for forced wins (Mate in 2) ---
        winning_moves_m2 = []
        for o_col in possible_moves:
            # Simulate O's move
            board_after_o = copy.deepcopy(self.board)
            o_row = self.get_next_open_row(board_after_o, o_col)
            board_after_o[o_row][o_col] = self.PLAYER_O

            # Assume X will make a move. O wins if for ALL of X's moves, O has a winning response.
            x_possible_moves = self.get_valid_locations(board_after_o)
            if not x_possible_moves: # No moves for X, so not a useful path
                continue
            
            o_wins_all_scenarios = True
            for x_col in x_possible_moves:
                # Simulate X's move
                board_after_x = copy.deepcopy(board_after_o)
                x_row = self.get_next_open_row(board_after_x, x_col)
                board_after_x[x_row][x_col] = self.PLAYER_X

                # Does O have a winning move now?
                o_can_win_now = False
                o_next_moves = self.get_valid_locations(board_after_x)
                for o_next_col in o_next_moves:
                    board_final = copy.deepcopy(board_after_x)
                    o_next_row = self.get_next_open_row(board_final, o_next_col)
                    board_final[o_next_row][o_next_col] = self.PLAYER_O
                    if self.check_win(board_final, self.PLAYER_O):
                        o_can_win_now = True
                        break # Found a winning move for O
                
                if not o_can_win_now:
                    o_wins_all_scenarios = False
                    break # X found a move that O cannot counter with a win

            if o_wins_all_scenarios:
                move_str = f"{chr(ord('a') + o_col)}{self.ROWS - o_row}"
                winning_moves_m2.append(move_str)

        if winning_moves_m2:
            print(", ".join(sorted(winning_moves_m2)))
            return
            
        print("No optimal winning moves found.")


# Board state from the problem:
# . . . . . . .
# . . . . . . .
# . . . . . . .
# . . . O O . .
# O . X O X X X
# X O O X X O X
#
# 1 represents 'X', 2 represents 'O', 0 is empty.
initial_board = [
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 2, 0, 0],
    [2, 0, 1, 2, 1, 1, 1],
    [1, 2, 2, 1, 1, 2, 1]
]

solver = Connect4Solver(initial_board)
solver.solve()