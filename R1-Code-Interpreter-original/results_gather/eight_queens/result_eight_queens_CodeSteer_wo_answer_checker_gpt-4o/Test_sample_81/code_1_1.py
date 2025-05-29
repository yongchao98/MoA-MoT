def is_safe(board, row, col, n):
    # Check if the current column has a queen
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper left diagonal
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper right diagonal
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False

    # Check if the current position is restricted
    if board[row][col] == 'X':
        return False

    return True

def solve_n_queens(board, row, n, result):
    if row >= n:
        return True

    for col in range(n):
        # Check if the position is safe and not restricted
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            result.append((row, col))
            if solve_n_queens(board, row + 1, n, result):
                return True
            # Backtrack
            board[row][col] = 0
            result.pop()

    return False

def place_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 'X', 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 'X', 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    result = []
    if solve_n_queens(board, 0, n, result):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in result) + ">>>")
    else:
        print("No solution found")

place_queens()