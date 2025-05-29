def is_safe(board, row, col):
    # Check this row on left side
    for i in range(col):
        if board[row][i] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check lower diagonal on left side
    for i, j in zip(range(row, len(board), 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens_util(board, row):
    # Base case: If all queens are placed
    if row >= len(board):
        return True

    # If there's already a queen in this row, move to the next row
    if 1 in board[row]:
        return solve_n_queens_util(board, row + 1)

    # Consider this row and try placing this queen in all columns one by one
    for col in range(len(board)):
        if board[row][col] == 'X' or board[row][col] == 1:  # Skip 'X' and already placed queens
            continue
        if is_safe(board, row, col):
            board[row][col] = 1
            if solve_n_queens_util(board, row + 1):
                return True
            board[row][col] = 0  # Backtrack

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

    # Convert 'X' to a special marker
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 'X':
                board[i][j] = 'X'
            elif board[i][j] == 1:
                board[i][j] = 1
            else:
                board[i][j] = 0

    if not solve_n_queens_util(board, 0):
        return "Solution does not exist"

    # Collect the positions of the queens
    queen_positions = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 1:
                queen_positions.append(f"{i} {j}")

    return queen_positions

# Print the solution
print(solve_n_queens())