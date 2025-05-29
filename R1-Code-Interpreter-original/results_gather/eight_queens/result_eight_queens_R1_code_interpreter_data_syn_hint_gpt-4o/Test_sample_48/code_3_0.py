def is_safe(board, row, col, n, cols, diag1, diag2):
    return not (cols[col] or diag1[row - col] or diag2[row + col])

def solve_n_queens(board, row, n, positions, cols, diag1, diag2):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(board, row, col, n, cols, diag1, diag2):
            board[row][col] = 1
            positions.append((row, col))
            cols[col] = diag1[row - col] = diag2[row + col] = True

            if solve_n_queens(board, row + 1, n, positions, cols, diag1, diag2):
                return True

            board[row][col] = 0
            positions.pop()
            cols[col] = diag1[row - col] = diag2[row + col] = False

    return False

def find_queen_positions():
    board = [
        ['X', 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 'X', 0, 0, 0, 0, 0],
        [0, 'X', 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    n = 8
    positions = [(1, 6), (3, 2)]
    for row, col in positions:
        board[row][col] = 1

    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)

    for row, col in positions:
        cols[col] = diag1[row - col] = diag2[row + col] = True

    start_row = 0
    while start_row < n and any(board[start_row][col] == 1 for col in range(n)):
        start_row += 1

    solve_n_queens(board, start_row, n, positions, cols, diag1, diag2)

    return positions

queen_positions = find_queen_positions()
print(queen_positions)