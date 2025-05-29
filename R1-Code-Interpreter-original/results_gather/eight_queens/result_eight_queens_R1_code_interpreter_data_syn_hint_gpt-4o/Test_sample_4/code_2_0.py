def is_safe(row, col, columns, diag1, diag2):
    return col not in columns and (row - col) not in diag1 and (row + col) not in diag2

def solve_n_queens_util(n, row, columns, diag1, diag2, board):
    if row == n:
        return True

    if row == 2:  # Skip row 2 as it already has a queen
        return solve_n_queens_util(n, row + 1, columns, diag1, diag2, board)

    for col in range(n):
        if is_safe(row, col, columns, diag1, diag2):
            board[row][col] = 1
            columns.add(col)
            diag1.add(row - col)
            diag2.add(row + col)

            if solve_n_queens_util(n, row + 1, columns, diag1, diag2, board):
                return True

            board[row][col] = 0
            columns.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)

    return False

def solve_n_queens():
    n = 8
    board = [[0 for _ in range(n)] for _ in range(n)]
    board[2][0] = 1  # Initial queen placement

    columns = {0}  # Column 0 is already occupied
    diag1 = {2 - 0}  # Diagonal occupied by the initial queen
    diag2 = {2 + 0}  # Diagonal occupied by the initial queen

    if not solve_n_queens_util(n, 0, columns, diag1, diag2, board):
        return "No solution exists"

    result = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                result.append(f"{i} {j}")
    return result

print(solve_n_queens())