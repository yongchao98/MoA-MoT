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
        if board[row][i] == 'X':
            continue
        if is_safe(board, row, i, n):
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

def solve_n_queens():
    n = 8
    board = [[0 for _ in range(n)] for _ in range(n)]
    board[3][3] = 'X'
    board[5][1] = 'X'

    if not solve_n_queens_util(board, 0, n):
        return "No solution exists"

    # Collect the positions of the queens
    queen_positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queen_positions.append(f"{i} {j}")

    return queen_positions

# Execute the function and print the result
print(solve_n_queens())