def is_safe(board, row, col):
    # Check this column on upper side
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper diagonal on right side
    for i, j in zip(range(row, -1, -1), range(col, len(board))):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row):
    if row == 7:  # Skip row 7 as it already has a queen
        return True

    for col in range(len(board)):
        if is_safe(board, row, col):
            board[row][col] = 1
            if solve_n_queens(board, row + 1):
                return True
            board[row][col] = 0

    return False

# Initialize the board
board = [[0 for _ in range(8)] for _ in range(8)]
board[7][7] = 1  # Place the initial queen

# Solve the problem
solve_n_queens(board, 0)

# Collect the positions of the queens
queen_positions = []
for i in range(8):
    for j in range(8):
        if board[i][j] == 1:
            queen_positions.append((i, j))

print(queen_positions)