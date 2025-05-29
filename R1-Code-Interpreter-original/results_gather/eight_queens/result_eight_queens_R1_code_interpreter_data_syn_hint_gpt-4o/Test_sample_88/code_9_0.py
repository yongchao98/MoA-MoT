def is_safe(board, row, col):
    # Check this column on upper side
    for i in range(row):
        if board[i] == col or \
           board[i] - i == col - row or \
           board[i] + i == col + row:
            return False
    return True

def solve_n_queens(board, row, fixed_positions):
    if row == len(board):
        return True
    if row in fixed_positions:
        return solve_n_queens(board, row + 1, fixed_positions)
    for col in range(len(board)):
        if is_safe(board, row, col):
            board[row] = col
            if solve_n_queens(board, row + 1, fixed_positions):
                return True
            board[row] = -1
    return False

# Initialize the board with -1 (no queen placed)
board = [-1] * 8

# Place the queens that are already on the board
fixed_positions = {3: 2, 6: 7}
for row, col in fixed_positions.items():
    board[row] = col

# Solve the problem
solve_n_queens(board, 0, fixed_positions)

# Collect the positions of the queens
queen_positions = [(i, board[i]) for i in range(8)]
print(queen_positions)