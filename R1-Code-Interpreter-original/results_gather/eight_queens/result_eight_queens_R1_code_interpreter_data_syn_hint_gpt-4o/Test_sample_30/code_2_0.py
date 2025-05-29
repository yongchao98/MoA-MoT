def is_safe(row, col, cols, diag1, diag2):
    return not (cols[col] or diag1[row - col] or diag2[row + col])

def place_queens(board, row, n, cols, diag1, diag2, positions):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 'X' or board[row][col] == 1:
            continue
        if is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            cols[col] = diag1[row - col] = diag2[row + col] = True
            positions.append((row, col))
            if place_queens(board, row + 1, n, cols, diag1, diag2, positions):
                return True
            board[row][col] = 0
            cols[col] = diag1[row - col] = diag2[row + col] = False
            positions.pop()

    return False

def find_queen_positions():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 'X', 0]
    ]

    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)

    # Mark the initial queens and restricted positions
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                cols[j] = diag1[i - j] = diag2[i + j] = True
                positions.append((i, j))

    place_queens(board, 0, n, cols, diag1, diag2, positions)

    return [f"{i} {j}" for i, j in positions]

print(find_queen_positions())