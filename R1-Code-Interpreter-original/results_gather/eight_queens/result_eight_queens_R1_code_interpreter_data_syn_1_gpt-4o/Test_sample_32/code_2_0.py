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

def solve_n_queens_util(board, row, n, positions):
    if row >= n:
        return True

    # If there's already a queen in this row, skip to the next row
    if any(board[row][i] == 1 for i in range(n)):
        return solve_n_queens_util(board, row + 1, n, positions)

    for i in range(n):
        if board[row][i] == 0 and is_safe(board, row, i, n):
            board[row][i] = 1
            positions.append((row, i))

            if solve_n_queens_util(board, row + 1, n, positions):
                return True

            board[row][i] = 0
            positions.pop()

    return False

def solve_n_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 'X', 0, 0, 0, 0, 0]
    ]

    positions = [(4, 4), (6, 0)]

    if solve_n_queens_util(board, 0, n, positions):
        return positions
    else:
        return "No solution exists"

positions = solve_n_queens()
print(positions)