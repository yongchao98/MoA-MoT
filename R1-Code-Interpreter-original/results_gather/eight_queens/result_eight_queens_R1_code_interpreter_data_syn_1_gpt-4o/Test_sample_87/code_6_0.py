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

    for col in range(n):
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            result.append((row, col))

            if solve_n_queens_util(board, row + 1, n, result):
                return True

            board[row][col] = 0
            result.pop()

    return False

def solve_n_queens(board):
    n = len(board)
    result = []

    # Mark initial queens
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                result.append((i, j))

    # Start solving from the first row that does not have a queen
    start_row = 0
    for i in range(n):
        if any(board[i][j] == 1 for j in range(n)):
            start_row = i + 1

    if not solve_n_queens_util(board, start_row, n, result):
        return "No solution exists"

    return result

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, -1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, -1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, -1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# Solve the problem
solution = solve_n_queens(board)

# Format the solution
formatted_solution = ', '.join(f"{r} {c}" for r, c in solution)
print(f"<<<{formatted_solution}>>>")