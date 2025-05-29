def solve_n_queens_optimized(board, row, n, cols, main_diag, anti_diag, result):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and col not in cols and (row - col) not in main_diag and (row + col) not in anti_diag:
            # Place the queen
            board[row][col] = 1
            cols.add(col)
            main_diag.add(row - col)
            anti_diag.add(row + col)
            result.append((row, col))

            # Recurse to place the next queen
            if solve_n_queens_optimized(board, row + 1, n, cols, main_diag, anti_diag, result):
                return True

            # Backtrack
            board[row][col] = 0
            cols.remove(col)
            main_diag.remove(row - col)
            anti_diag.remove(row + col)
            result.pop()

    return False

def place_queens_optimized():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 'X', 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 'X', 0, 0, 0, 'X', 'X', 0]
    ]

    # Convert 'X' to -1 for easier processing
    for i in range(n):
        for j in range(n):
            if board[i][j] == 'X':
                board[i][j] = -1

    cols = set()
    main_diag = set()
    anti_diag = set()
    result = []

    if solve_n_queens_optimized(board, 0, n, cols, main_diag, anti_diag, result):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in result) + ">>>")
    else:
        print("No solution found")

place_queens_optimized()