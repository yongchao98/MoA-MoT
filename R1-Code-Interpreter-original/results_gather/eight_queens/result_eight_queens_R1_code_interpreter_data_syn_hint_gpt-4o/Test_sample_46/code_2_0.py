def is_safe(row, col, cols, diag1, diag2, n):
    return not (cols[col] or diag1[row - col + n - 1] or diag2[row + col])

def solve_n_queens(board, row, n, cols, diag1, diag2):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(row, col, cols, diag1, diag2, n):
            board[row][col] = 1
            cols[col] = diag1[row - col + n - 1] = diag2[row + col] = True
            if solve_n_queens(board, row + 1, n, cols, diag1, diag2):
                return True
            board[row][col] = 0
            cols[col] = diag1[row - col + n - 1] = diag2[row + col] = False

    return False

def find_queen_positions(board):
    n = len(board)
    cols = [False] * n
    diag1 = [False] * (2 * n - 1)
    diag2 = [False] * (2 * n - 1)

    # Mark the already placed queens and restricted positions
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                cols[j] = diag1[i - j + n - 1] = diag2[i + j] = True
            elif board[i][j] == 'X':
                board[i][j] = -1

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
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 'X'],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 'X', 0, 0, 0, 0, 0]
]

# Find and print the positions of the queens
queen_positions = find_queen_positions(board)
print("<<<" + ", ".join(queen_positions) + ">>>")