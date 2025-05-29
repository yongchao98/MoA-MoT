def is_safe(row, col, cols, diag1, diag2, board):
    # Check if the position is not blocked by 'X' and is not under attack
    return (col not in cols and
            (row - col) not in diag1 and
            (row + col) not in diag2 and
            board[row][col] != -1)

def solve_n_queens_optimized(board, row, n, positions, cols, diag1, diag2):
    if row >= n:
        return True

    for col in range(n):
        if is_safe(row, col, cols, diag1, diag2, board):
            # Place the queen
            board[row][col] = 1
            positions.append((row, col))
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)

            if solve_n_queens_optimized(board, row + 1, n, positions, cols, diag1, diag2):
                return True

            # Backtrack
            board[row][col] = 0
            positions.pop()
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)

    return False

def place_queens_optimized():
    board = [
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 'X', 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    n = 8
    positions = []

    # Convert 'X' to a non-placeable marker
    for i in range(n):
        for j in range(n):
            if board[i][j] == 'X':
                board[i][j] = -1

    # Sets to track columns and diagonals
    cols = set()
    diag1 = set()
    diag2 = set()

    # Start solving from the first row
    if solve_n_queens_optimized(board, 0, n, positions, cols, diag1, diag2):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
    else:
        print("No solution found")

place_queens_optimized()