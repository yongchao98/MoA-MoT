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

    # Skip row 2 as it already has a queen
    if row == 2:
        return solve_n_queens_util(board, row + 1, n)

    for i in range(n):
        if is_safe(board, row, i, n):
            board[row][i] = 1
            if solve_n_queens_util(board, row + 1, n):
                return True
            board[row][i] = 0

    return False

def solve_n_queens():
    n = 8
    board = [[0 for _ in range(n)] for _ in range(n)]
    board[2][0] = 1  # Initial queen placement

    if not solve_n_queens_util(board, 0, n):
        return "No solution exists"

    result = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                result.append(f"{i} {j}")
    return result

print(solve_n_queens())