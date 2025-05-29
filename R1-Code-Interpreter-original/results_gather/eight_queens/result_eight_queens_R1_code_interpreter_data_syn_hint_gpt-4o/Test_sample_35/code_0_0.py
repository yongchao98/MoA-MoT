def is_safe(board, row, col, n):
    # Check this column on upper side
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper diagonal on right side
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, n, positions):
    if row >= n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            positions.append((row, col))
            if solve_n_queens(board, row + 1, n, positions):
                return True
            board[row][col] = 0
            positions.pop()

    return False

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 'X', 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 'X'],
    [0, 0, 0, 'X', 0, 0, 0, 0]
]

# Convert 'X' to -1 for easier processing
for i in range(8):
    for j in range(8):
        if board[i][j] == 'X':
            board[i][j] = -1

positions = [(1, 0), (2, 3)]
solve_n_queens(board, 0, 8, positions)

# Format the output
output = ', '.join(f"{r} {c}" for r, c in positions)
print(f"<<<{output}>>>")