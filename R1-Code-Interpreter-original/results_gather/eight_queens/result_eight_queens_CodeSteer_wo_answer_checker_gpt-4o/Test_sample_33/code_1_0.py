def is_safe(row, col, cols, diag1, diag2, board):
    return not (cols[col] or diag1[row - col + len(board) - 1] or diag2[row + col] or board[row][col] == 'X')

def solve_n_queens(board, row, n, result, cols, diag1, diag2):
    if row == n:
        return True

    for col in range(n):
        if is_safe(row, col, cols, diag1, diag2, board):
            board[row][col] = 1
            result.append((row, col))
            cols[col] = diag1[row - col + n - 1] = diag2[row + col] = True

            if solve_n_queens(board, row + 1, n, result, cols, diag1, diag2):
                return True

            board[row][col] = 0
            result.pop()
            cols[col] = diag1[row - col + n - 1] = diag2[row + col] = False

    return False

def place_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    result = []
    cols = [False] * n
    diag1 = [False] * (2 * n - 1)
    diag2 = [False] * (2 * n - 1)

    # Mark initial queens
    initial_queens = [(0, 4), (4, 3)]
    for r, c in initial_queens:
        board[r][c] = 1
        cols[c] = diag1[r - c + n - 1] = diag2[r + c] = True
        result.append((r, c))

    # Start placing queens from the first row that does not have a queen
    start_row = 0
    while start_row < n and any(board[start_row]):
        start_row += 1

    solve_n_queens(board, start_row, n, result, cols, diag1, diag2)

    # Format the result
    formatted_result = ', '.join(f'{r} {c}' for r, c in result)
    print(f'<<<{formatted_result}>>>')

place_queens()