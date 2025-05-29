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
    if row >= len(board):
        return True

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
    [0, 0, 'X', 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 'X'],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# Convert 'X' to -1 for easier processing
for i in range(len(board)):
    for j in range(len(board[i])):
        if board[i][j] == 'X':
            board[i][j] = -1

# Initial positions of queens
queens_positions = [(1, 1), (6, 4)]

# Solve the problem
solve_n_queens(board, 0, queens_positions)

# Format the output
output = ', '.join(f"{r} {c}" for r, c in sorted(queens_positions))
print(f"<<<{output}>>>")