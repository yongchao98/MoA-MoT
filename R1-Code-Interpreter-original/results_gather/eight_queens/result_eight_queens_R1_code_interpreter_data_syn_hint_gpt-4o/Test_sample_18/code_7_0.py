def is_safe(row, col, cols, diag1, diag2, board):
    if col in cols or (row - col) in diag1 or (row + col) in diag2 or board[row][col] == 'X':
        return False
    return True

def solve_n_queens_util(row, n, board, cols, diag1, diag2, result):
    if row == n:
        return True

    for col in range(n):
        if is_safe(row, col, cols, diag1, diag2, board):
            board[row][col] = 1
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)
            result.append((row, col))

            if solve_n_queens_util(row + 1, n, board, cols, diag1, diag2, result):
                return True

            board[row][col] = 0
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)
            result.pop()

    return False

def solve_n_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 0, 0, 'X', 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        ['X', 0, 0, 0, 0, 0, 0, 0]
    ]

    cols = set()
    diag1 = set()
    diag2 = set()
    result = []

    if solve_n_queens_util(0, n, board, cols, diag1, diag2, result):
        return result
    else:
        return "No solution exists"

solution = solve_n_queens()
print(solution)