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
    # Base case: If all queens are placed
    if row >= n:
        return True

    # Consider this row and try placing this queen in all columns one by one
    for i in range(n):
        if board[row][i] == 0 and is_safe(board, row, i, n):
            # Place this queen in board[row][i]
            board[row][i] = 1

            # Recur to place rest of the queens
            if solve_n_queens_util(board, row + 1, n):
                return True

            # If placing queen in board[row][i] doesn't lead to a solution
            # then remove queen from board[row][i]
            board[row][i] = 0

    # If the queen cannot be placed in any column in this row
    return False

def solve_n_queens(n):
    board = [[0 for _ in range(n)] for _ in range(n)]

    if not solve_n_queens_util(board, 0, n):
        return "Solution does not exist"

    # Collect the positions of the queens
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Define the board size
n = 8

# Solve the N-Queens problem
queen_positions = solve_n_queens(n)
print(queen_positions)