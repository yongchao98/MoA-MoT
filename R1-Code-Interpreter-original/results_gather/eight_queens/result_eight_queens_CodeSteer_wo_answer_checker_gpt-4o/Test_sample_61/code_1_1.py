def solve_n_queens_optimized(board, row, n, cols, diag1, diag2, result):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == -1:  # Skip 'X' positions
            continue
        if col in cols or (row - col) in diag1 or (row + col) in diag2:
            continue

        # Place the queen
        board[row][col] = 1
        cols.add(col)
        diag1.add(row - col)
        diag2.add(row + col)
        result.append((row, col))

        if solve_n_queens_optimized(board, row + 1, n, cols, diag1, diag2, result):
            return True

        # Backtrack
        board[row][col] = 0
        cols.remove(col)
        diag1.remove(row - col)
        diag2.remove(row + col)
        result.pop()

    return False

def place_queens_optimized():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, -1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    n = 8
    result = []
    cols = set()
    diag1 = set()
    diag2 = set()

    if solve_n_queens_optimized(board, 0, n, cols, diag1, diag2, result):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in result) + ">>>")
    else:
        print("No solution found")

place_queens_optimized()