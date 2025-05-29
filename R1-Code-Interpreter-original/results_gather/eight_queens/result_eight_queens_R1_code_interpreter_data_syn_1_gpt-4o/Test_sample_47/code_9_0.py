def is_safe(row, col, occupied_rows, occupied_cols, occupied_diag1, occupied_diag2):
    if row in occupied_rows or col in occupied_cols or (row - col) in occupied_diag1 or (row + col) in occupied_diag2:
        return False
    return True

def solve_n_queens_util(board, col, n, result, occupied_rows, occupied_cols, occupied_diag1, occupied_diag2):
    # base case: If all queens are placed
    if col >= n:
        return True

    # If this column already has a queen, move to the next column
    for i in range(n):
        if board[i][col] == 1:
            return solve_n_queens_util(board, col + 1, n, result, occupied_rows, occupied_cols, occupied_diag1, occupied_diag2)

    # Consider this column and try placing this queen in all rows one by one
    for row in range(n):
        if board[row][col] == 'X' or board[row][col] == 1:
            continue
        if is_safe(row, col, occupied_rows, occupied_cols, occupied_diag1, occupied_diag2):
            # Place this queen in board[row][col]
            board[row][col] = 1
            result.append((row, col))
            occupied_rows.add(row)
            occupied_cols.add(col)
            occupied_diag1.add(row - col)
            occupied_diag2.add(row + col)

            # Recur to place rest of the queens
            if solve_n_queens_util(board, col + 1, n, result, occupied_rows, occupied_cols, occupied_diag1, occupied_diag2):
                return True

            # If placing queen in board[row][col] doesn't lead to a solution
            # then remove queen from board[row][col]
            board[row][col] = 0
            result.pop()
            occupied_rows.remove(row)
            occupied_cols.remove(col)
            occupied_diag1.remove(row - col)
            occupied_diag2.remove(row + col)

    # If the queen cannot be placed in any row in this column col then return false
    return False

def solve_n_queens(board):
    n = len(board)
    result = []
    occupied_rows = set()
    occupied_cols = set()
    occupied_diag1 = set()
    occupied_diag2 = set()

    # Mark existing queens in the result and occupied sets
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                result.append((i, j))
                occupied_rows.add(i)
                occupied_cols.add(j)
                occupied_diag1.add(i - j)
                occupied_diag2.add(i + j)

    # Start from the first column
    if not solve_n_queens_util(board, 0, n, result, occupied_rows, occupied_cols, occupied_diag1, occupied_diag2):
        return "Solution does not exist"
    return result

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 'X'],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    ['X', 0, 0, 'X', 0, 0, 0, 0],
    [0, 0, 0, 'X', 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# Solve the problem
solution = solve_n_queens(board)
print(solution)