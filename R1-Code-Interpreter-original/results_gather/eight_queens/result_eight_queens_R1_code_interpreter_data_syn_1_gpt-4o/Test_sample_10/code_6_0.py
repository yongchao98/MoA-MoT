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

def solve_n_queens_util(board, col, queens_placed):
    # Base case: If all queens are placed
    if queens_placed == 8:
        return True

    # Consider this column and try placing this queen in all rows one by one
    for i in range(len(board)):
        if board[i][col] == 'X' or board[i][col] == 1:  # Skip the blocked cell and already placed queen
            continue
        if is_safe(board, i, col):
            # Place this queen in board[i][col]
            board[i][col] = 1

            # Recur to place rest of the queens
            if solve_n_queens_util(board, col + 1, queens_placed + 1):
                return True

            # If placing queen in board[i][col] doesn't lead to a solution
            # then remove queen from board[i][col]
            board[i][col] = 0

    # If the queen cannot be placed in any row in this column col then return false
    return False

def solve_n_queens():
    # Initialize the board
    board = [[0, 0, 0, 0, 0, 0, 0, 0] for _ in range(8)]
    board[4][5] = 1  # Place the initial queen
    board[5][5] = 'X'  # Mark the blocked cell

    # Start placing queens from the first column
    if not solve_n_queens_util(board, 0, 1):  # 1 queen is already placed
        return "Solution does not exist"

    # Collect the positions of the queens
    queen_positions = []
    for i in range(8):
        for j in range(8):
            if board[i][j] == 1:
                queen_positions.append(f"{i} {j}")

    return queen_positions

# Get the solution
solution = solve_n_queens()
print(solution)