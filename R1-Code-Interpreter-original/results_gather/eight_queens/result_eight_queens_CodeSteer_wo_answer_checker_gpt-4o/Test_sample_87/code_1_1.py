def is_safe(row, col, cols, diag1, diag2):
    return col not in cols and (row - col) not in diag1 and (row + col) not in diag2

def solve_n_queens(board, row, n, positions, cols, diag1, diag2):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            positions.append((row, col))
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)

            if solve_n_queens(board, row + 1, n, positions, cols, diag1, diag2):
                return True

            # Backtrack
            board[row][col] = 0
            positions.pop()
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)

    return False

def place_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 'X', 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    positions = [(5, 4)]
    cols = {4}
    diag1 = {5 - 4}
    diag2 = {5 + 4}

    if solve_n_queens(board, 0, n, positions, cols, diag1, diag2):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
    else:
        print("No solution found")

place_queens()