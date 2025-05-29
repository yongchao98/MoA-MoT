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

def solve_n_queens(board, row, n):
    if row >= n:
        return True

    # If the row already has a queen, move to the next row
    if any(board[row][j] == 1 for j in range(n)):
        return solve_n_queens(board, row + 1, n)

    for col in range(n):
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            if solve_n_queens(board, row + 1, n):
                return True
            board[row][col] = 0

    return False

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 1, -1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [-1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [-1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 0, 0, 0, 0]
]

n = len(board)

# Solve the board
solve_n_queens(board, 0, n)

# Collect the positions of the queens
queen_positions = []
for i in range(n):
    for j in range(n):
        if board[i][j] == 1:
            queen_positions.append((i, j))

# Format the output
output = ', '.join(f"{r} {c}" for r, c in queen_positions)
print(f"<<<{output}>>>")