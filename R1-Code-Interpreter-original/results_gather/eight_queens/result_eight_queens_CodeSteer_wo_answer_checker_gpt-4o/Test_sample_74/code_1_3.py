def is_safe(row, col, cols, diag1, diag2):
    return col not in cols and (row - col) not in diag1 and (row + col) not in diag2

def solve_n_queens(board, row, n, cols, diag1, diag2, positions):
    if row == n:
        return True

    # Skip rows that already have queens
    if any(board[row]):
        return solve_n_queens(board, row + 1, n, cols, diag1, diag2, positions)

    for col in range(n):
        if board[row][col] == 0 and is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)
            positions.append((row, col))

            if solve_n_queens(board, row + 1, n, cols, diag1, diag2, positions):
                return True

            board[row][col] = 0
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)
            positions.pop()

    return False

def place_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0]
    ]

    # Initial positions of queens
    positions = [(1, 4), (7, 0)]
    # Track columns and diagonals occupied by queens
    cols = {4, 0}
    diag1 = {1 - 4, 7 - 0}
    diag2 = {1 + 4, 7 + 0}

    # Start solving from the first row
    if solve_n_queens(board, 0, n, cols, diag1, diag2, positions):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
    else:
        print("No solution found")

place_queens()