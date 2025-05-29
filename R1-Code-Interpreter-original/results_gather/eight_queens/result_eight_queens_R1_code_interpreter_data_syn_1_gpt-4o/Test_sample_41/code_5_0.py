def is_safe(row, col, rows, cols, diag1, diag2, n):
    return not (rows[row] or cols[col] or diag1[row - col + n - 1] or diag2[row + col])

def place_queens(board, row, rows, cols, diag1, diag2, n, positions):
    if len(positions) == n:
        return True

    for col in range(n):
        if board[row][col] == 'X' or board[row][col] == 1:
            continue
        if is_safe(row, col, rows, cols, diag1, diag2, n):
            board[row][col] = 1
            rows[row] = cols[col] = diag1[row - col + n - 1] = diag2[row + col] = True
            positions.append((row, col))

            if place_queens(board, row + 1, rows, cols, diag1, diag2, n, positions):
                return True

            board[row][col] = 0
            rows[row] = cols[col] = diag1[row - col + n - 1] = diag2[row + col] = False
            positions.pop()

    return False

def solve_n_queens(board):
    n = len(board)
    rows = [False] * n
    cols = [False] * n
    diag1 = [False] * (2 * n - 1)
    diag2 = [False] * (2 * n - 1)
    positions = []

    # Mark the initial positions of the queens
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                rows[i] = cols[j] = diag1[i - j + n - 1] = diag2[i + j] = True
                positions.append((i, j))

    # Start placing queens from the first unoccupied row
    start_row = 0
    for i in range(n):
        if not rows[i]:
            start_row = i
            break

    if not place_queens(board, start_row, rows, cols, diag1, diag2, n, positions):
        return "Solution does not exist"
    else:
        return [f"{r} {c}" for r, c in positions]

# Initial board setup
board = [
    ['X', 'X', 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 'X'],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 'X', 0, 0, 0, 0, 0]
]

# Solve the N-Queens problem
queen_positions = solve_n_queens(board)
print(queen_positions)