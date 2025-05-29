def is_valid(board, row, col, n):
    # Check if the column is safe
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check the upper left diagonal
    i, j = row, col
    while i >= 0 and j >= 0:
        if board[i][j] == 1:
            return False
        i -= 1
        j -= 1

    # Check the upper right diagonal
    i, j = row, col
    while i >= 0 and j < n:
        if board[i][j] == 1:
            return False
        i -= 1
        j += 1

    return True

def solve_n_queens(board, row, n, positions):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_valid(board, row, col, n):
            board[row][col] = 1
            positions.append((row, col))
            if solve_n_queens(board, row + 1, n, positions):
                return True
            board[row][col] = 0
            positions.pop()

    return False

def place_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 'X', 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 'X', 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    positions = []
    if solve_n_queens(board, 0, n, positions):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
    else:
        print("No solution found")

place_queens()