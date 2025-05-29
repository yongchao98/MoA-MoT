def is_safe(board, row, col):
    # Check this column on upper side
    for i in range(row):
        if board[i] == col or \
           board[i] - i == col - row or \
           board[i] + i == col + row:
            return False
    return True

def solve_n_queens(board, row):
    if row == len(board):
        return True
    for col in range(len(board)):
        if is_safe(board, row, col):
            board[row] = col
            if solve_n_queens(board, row + 1):
                return True
            board[row] = -1
    return False

# Initialize the board with -1 (no queen placed)
board = [-1] * 8

# Place the queens that are already on the board
board[3] = 2
board[6] = 7

# Solve the problem
solve_n_queens(board, 0)

# Collect the positions of the queens
queen_positions = [(i, board[i]) for i in range(8)]
print(queen_positions)