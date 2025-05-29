def is_safe(board, row, col):
    # Check this column on upper side
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper diagonal on right side
    for i, j in zip(range(row, -1, -1), range(col, len(board))):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, queens_positions):
    if row >= len(board):
        return True

    for col in range(len(board)):
        if is_safe(board, row, col):
            board[row][col] = 1
            queens_positions.append((row, col))

            if solve_n_queens(board, row + 1, queens_positions):
                return True

            board[row][col] = 0
            queens_positions.pop()

    return False

# Initialize the board
board = [[0] * 8 for _ in range(8)]
board[3][3] = 1  # Place the initial queen

# List to store the positions of queens
queens_positions = [(3, 3)]

# Solve the problem
solve_n_queens(board, 0, queens_positions)

# Print the positions of the queens
print(queens_positions)