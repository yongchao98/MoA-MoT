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

def solve_n_queens_util(board, row, n):
    if row >= n:
        return True

    for i in range(n):
        if board[row][i] == 0 and is_safe(board, row, i, n):
            board[row][i] = 1
            if solve_n_queens_util(board, row + 1, n):
                return True
            board[row][i] = 0

    return False

def solve_n_queens(n):
    board = [[0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 'X', 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 'X'],
             [0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0]]

    for i in range(n):
        for j in range(n):
            if board[i][j] == 'X':
                board[i][j] = -1

    if not solve_n_queens_util(board, 0, n):
        return "No solution exists"

    result = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                result.append(f"{i} {j}")

    return result

print(solve_n_queens(8))