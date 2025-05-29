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

def solve_n_queens(board, row, queens_positions):
    if len(queens_positions) == 8:
        return True

    if row == 4:  # Skip the row where a queen is already placed
        return solve_n_queens(board, row + 1, queens_positions)

    for col in range(len(board)):
        if board[row][col] == 0 and is_safe(board, row, col):
            board[row][col] = 1
            queens_positions.append((row, col))

            if solve_n_queens(board, row + 1, queens_positions):
                return True

            board[row][col] = 0
            queens_positions.pop()

    return False

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

queens_positions = [(4, 2)]

# Solve the problem starting from the first row
solve_n_queens(board, 0, queens_positions)

# Output the positions of the queens
print(queens_positions)