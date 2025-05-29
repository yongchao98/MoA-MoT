def is_safe(row, col, cols, diag1, diag2):
    return col not in cols and (row - col) not in diag1 and (row + col) not in diag2

def solve_n_queens_optimized(board, row, queens, cols, diag1, diag2):
    if row == len(board):
        return True

    for col in range(len(board)):
        if board[row][col] == 'X':
            continue
        if is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            queens.append((row, col))
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)

            if solve_n_queens_optimized(board, row + 1, queens, cols, diag1, diag2):
                return True

            board[row][col] = 0
            queens.pop()
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)

    return False

def place_queens_optimized():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 'X', 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 'X', 0, 0, 0, 0, 0, 0]
    ]

    queens = []
    cols = set()
    diag1 = set()
    diag2 = set()

    solve_n_queens_optimized(board, 0, queens, cols, diag1, diag2)
    print("<<<" + ", ".join(f"{r} {c}" for r, c in queens) + ">>>")

place_queens_optimized()