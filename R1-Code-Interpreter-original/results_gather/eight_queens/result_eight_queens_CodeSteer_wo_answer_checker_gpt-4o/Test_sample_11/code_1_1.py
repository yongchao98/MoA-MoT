def is_safe(board, row, col, n):
    # Check if the position is blocked
    if board[row][col] == 'X':
        return False

    # Check this column on upper side
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper diagonal on right side
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, n, positions):
    if row >= n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            positions.append((row, col))

            if solve_n_queens(board, row + 1, n, positions):
                return True

            # Backtrack
            board[row][col] = 0
            positions.pop()

    return False

def place_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    positions = []
    if solve_n_queens(board, 0, n, positions):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
    else:
        print("No solution found")

place_queens()