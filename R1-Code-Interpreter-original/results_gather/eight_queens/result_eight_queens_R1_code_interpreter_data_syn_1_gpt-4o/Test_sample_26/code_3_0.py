def is_safe(row, col, columns, diag1, diag2):
    if columns[col] or diag1[row - col] or diag2[row + col]:
        return False
    return True

def solve_n_queens_util(board, row, columns, diag1, diag2):
    # Base case: If all queens are placed
    if row >= len(board):
        return True

    # If there's already a queen in this row, move to the next row
    if 1 in board[row]:
        return solve_n_queens_util(board, row + 1, columns, diag1, diag2)

    # Consider this row and try placing this queen in all columns one by one
    for col in range(len(board)):
        if board[row][col] == 'X' or board[row][col] == 1:  # Skip 'X' and already placed queens
            continue
        if is_safe(row, col, columns, diag1, diag2):
            board[row][col] = 1
            columns[col] = True
            diag1[row - col] = True
            diag2[row + col] = True
            if solve_n_queens_util(board, row + 1, columns, diag1, diag2):
                return True
            board[row][col] = 0  # Backtrack
            columns[col] = False
            diag1[row - col] = False
            diag2[row + col] = False

    return False

def solve_n_queens():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0]
    ]

    n = len(board)
    columns = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)

    # Mark the already placed queens
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                columns[j] = True
                diag1[i - j] = True
                diag2[i + j] = True

    if not solve_n_queens_util(board, 0, columns, diag1, diag2):
        return "Solution does not exist"

    # Collect the positions of the queens
    queen_positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queen_positions.append(f"{i} {j}")

    return queen_positions

# Print the solution
print(solve_n_queens())