def is_valid(board, row, col):
    # Check column
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, len(board))):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, positions):
    if row == len(board):
        return positions

    for col in range(len(board)):
        if board[row][col] == 0 and is_valid(board, row, col):
            board[row][col] = 1
            positions.append((row, col))
            result = solve_n_queens(board, row + 1, positions)
            if result:
                return result
            board[row][col] = 0
            positions.pop()

    return None

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [-1, 0, 0, 0, -1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0, 0, 0]
]

# Solve the problem
positions = solve_n_queens(board, 0, [])

# Format the output
if positions:
    output = ', '.join(f"{r} {c}" for r, c in positions)
    print(f"<<<{output}>>>")
else:
    print("No solution found.")