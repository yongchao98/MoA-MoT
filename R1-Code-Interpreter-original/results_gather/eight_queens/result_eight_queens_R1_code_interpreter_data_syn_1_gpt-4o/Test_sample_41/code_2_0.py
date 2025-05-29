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
    # base case: If all queens are placed
    if col >= n:
        return True

    # Consider this column and try placing this queen in all rows one by one
    for i in range(n):
        if board[i][col] == 'X' or board[i][col] == 1:  # Skip blocked cells and already placed queens
            continue
        if is_safe(board, i, col, n):
            # Place this queen in board[i][col]
            board[i][col] = 1

            # Recur to place rest of the queens
            if solve_n_queens_util(board, col + 1, n):
                return True

            # If placing queen in board[i][col] doesn't lead to a solution
            # then remove queen from board[i][col]
            board[i][col] = 0

    # If the queen cannot be placed in any row in this column col then return false
    return False

def solve_n_queens():
    n = 8
    board = [
        ['X', 'X', 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 'X'],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 'X', 0, 0, 0, 0, 0]
    ]

    # Start solving from the first column
    if not solve_n_queens_util(board, 0, n):
        return "Solution does not exist"

    # Collect the positions of the queens
    queen_positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                queen_positions.append(f"{i} {j}")

    return queen_positions

# Get the solution
solution = solve_n_queens()
print(f"<<<{', '.join(solution)}>>>")