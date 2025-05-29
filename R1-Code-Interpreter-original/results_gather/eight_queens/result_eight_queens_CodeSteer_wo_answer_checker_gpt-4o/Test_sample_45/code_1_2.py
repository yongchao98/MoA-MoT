def is_safe(board, row, col, n):
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

def solve_n_queens(board, row, n, result):
    if row >= n:
        return True

    for col in range(n):
        if board[row][col] == 'X' or board[row][col] == 1 or not is_safe(board, row, col, n):
            continue

        board[row][col] = 1
        result.append((row, col))
        print(f"Placing queen at: {row} {col}")

        if solve_n_queens(board, row + 1, n, result):
            return True

        # Backtrack
        board[row][col] = 0
        result.pop()
        print(f"Backtracking from: {row} {col}")

    return False

def place_queens():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 'X', 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    # Mark existing queens
    board[3][7] = 1

    result = []
    if solve_n_queens(board, 0, 8, result):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in result) + ">>>")
    else:
        print("No solution found")

place_queens()