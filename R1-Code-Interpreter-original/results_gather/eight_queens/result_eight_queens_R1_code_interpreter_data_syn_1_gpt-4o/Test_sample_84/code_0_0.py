def is_valid(board, row, col):
    # Check if there's a queen in the same column
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper left diagonal
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper right diagonal
    for i, j in zip(range(row, -1, -1), range(col, len(board))):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row, queens_positions):
    if row == len(board):
        return True

    for col in range(len(board)):
        if board[row][col] == 0 and is_valid(board, row, col):
            board[row][col] = 1
            queens_positions.append((row, col))
            if solve_n_queens(board, row + 1, queens_positions):
                return True
            board[row][col] = 0
            queens_positions.pop()

    return False

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    ['X', 0, 0, 0, 'X', 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 'X', 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 'X', 0, 0, 0, 0, 0]
]

# Convert 'X' to a non-placeable marker
for i in range(len(board)):
    for j in range(len(board[i])):
        if board[i][j] == 'X':
            board[i][j] = -1

queens_positions = []
solve_n_queens(board, 0, queens_positions)

# Convert positions to the required format
result = ', '.join(f"{r} {c}" for r, c in queens_positions)
print(f"<<<{result}>>>")