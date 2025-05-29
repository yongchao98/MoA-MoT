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

def solve_n_queens_util(board, col, n):
    # Base case: If all queens are placed
    if col >= n:
        return True

    # If there's already a queen in this column, move to the next column
    for i in range(n):
        if board[i][col] == 1:
            return solve_n_queens_util(board, col + 1, n)

    # Consider this column and try placing this queen in all rows one by one
    for i in range(n):
        if board[i][col] == 'X' or board[i][col] == 1:
            continue
        if is_safe(board, i, col, n):
            board[i][col] = 1
            if solve_n_queens_util(board, col + 1, n):
                return True
            board[i][col] = 0  # Backtrack

    return False

def solve_n_queens():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 'X', 0]
    ]

    # Start placing queens from the first column
    if not solve_n_queens_util(board, 0, n):
        return "Solution does not exist"

    # Collect the positions of the queens
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")

    return positions

# Print the solution
print(solve_n_queens())