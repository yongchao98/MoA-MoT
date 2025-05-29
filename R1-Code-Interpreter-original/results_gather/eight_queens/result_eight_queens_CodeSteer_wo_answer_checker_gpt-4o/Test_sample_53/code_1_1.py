def is_valid(board, row, col, n, x_positions):
    # Check if the position is marked with 'X'
    if (row, col) in x_positions:
        return False

    # Check column
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, n, positions, x_positions):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_valid(board, row, col, n, x_positions):
            board[row][col] = 1
            positions.append((row, col))
            if solve_n_queens(board, row + 1, n, positions, x_positions):
                return True
            board[row][col] = 0
            positions.pop()

    return False

def place_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    # Example 'X' positions
    x_positions = [(1, 2), (3, 5)]  # Add any 'X' positions here

    positions = [(0, 4)]
    solve_n_queens(board, 1, n, positions, x_positions)
    print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")

place_queens()