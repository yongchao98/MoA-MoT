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

def solve_n_queens_util(board, row, n, result):
    if row >= n:
        return True

    for i in range(n):
        if board[row][i] == 0 and is_safe(board, row, i, n):
            board[row][i] = 1
            result.append((row, i))

            if solve_n_queens_util(board, row + 1, n, result):
                return True

            board[row][i] = 0
            result.pop()

    return False

def solve_n_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 'X', 0, 0, 0, 0, 0],
        [0, 0, 'X', 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    # Convert 'X' to -1 for easier processing
    for i in range(n):
        for j in range(n):
            if board[i][j] == 'X':
                board[i][j] = -1

    result = []
    # Add already placed queens to result
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                result.append((i, j))

    if solve_n_queens_util(board, len(result), n, result):
        return result
    else:
        return "No solution exists"

solution = solve_n_queens()
print(solution)