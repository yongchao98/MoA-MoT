def is_safe(board, row, col, n, cols, diag1, diag2):
    if cols[col] or diag1[row - col] or diag2[row + col]:
        return False
    return True

def solve_n_queens(board, row, n, cols, diag1, diag2):
    if row >= n:
        return True

    for col in range(n):
        if board[row][col] == 'X' or board[row][col] == 1:
            continue
        if is_safe(board, row, col, n, cols, diag1, diag2):
            board[row][col] = 1
            cols[col] = diag1[row - col] = diag2[row + col] = True
            if solve_n_queens(board, row + 1, n, cols, diag1, diag2):
                return True
            board[row][col] = 0
            cols[col] = diag1[row - col] = diag2[row + col] = False

    return False

def find_queen_positions(board):
    n = len(board)
    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)

    # Mark the initial queens and 'X' positions
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                cols[j] = True
                diag1[i - j] = True
                diag2[i + j] = True

    # Start solving from the first row
    solve_n_queens(board, 0, n, cols, diag1, diag2)

    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    ['X', 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

positions = find_queen_positions(board)
print(positions)