def is_safe(board, row, col, n):
    # Check this row on left side
    for i in range(col):
        if board[row][i] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check lower diagonal on left side
    for i, j in zip(range(row, n, 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens_util(board, col, n, result):
    # base case: If all queens are placed
    if col >= n:
        return True

    # If there's already a queen in this column, skip it
    for i in range(n):
        if board[i][col] == 1:
            return solve_n_queens_util(board, col + 1, n, result)

    # Consider this column and try placing this queen in all rows one by one
    for i in range(n):
        if board[i][col] == 0 and is_safe(board, i, col, n):
            # Place this queen in board[i][col]
            board[i][col] = 1
            result.append((i, col))

            # Recur to place rest of the queens
            if solve_n_queens_util(board, col + 1, n, result):
                return True

            # If placing queen in board[i][col] doesn't lead to a solution
            # then remove queen from board[i][col]
            board[i][col] = 0
            result.pop()

    # If the queen cannot be placed in any row in this column col then return false
    return False

def solve_n_queens(board):
    n = len(board)
    result = []

    # Start from the first column and try to place queens
    if not solve_n_queens_util(board, 0, n, result):
        return "Solution does not exist"
    else:
        return result

# Initial board setup
board = [
    [0, 0, 0, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, -1],
    [0, -1, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# Solve the N-Queens problem
solution = solve_n_queens(board)

# Check if a solution was found
if isinstance(solution, str):
    print(solution)
else:
    # Format the solution
    formatted_solution = ', '.join(f"{r} {c}" for r, c in solution)
    print(f"<<<{formatted_solution}>>>")