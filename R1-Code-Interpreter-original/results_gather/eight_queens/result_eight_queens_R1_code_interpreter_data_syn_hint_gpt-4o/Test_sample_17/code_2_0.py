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

def solve_n_queens(board, row, queens_placed):
    if queens_placed == 8:
        return True

    if row == 3:  # Skip the row with the pre-placed queen
        return solve_n_queens(board, row + 1, queens_placed)

    for i in range(len(board)):
        if is_safe(board, row, i):
            board[row][i] = 1
            if solve_n_queens(board, row + 1, queens_placed + 1):
                return True
            board[row][i] = 0

    return False

def find_queen_positions(board):
    # Start solving from the first row
    queens_placed = sum(row.count(1) for row in board)
    solve_n_queens(board, 0, queens_placed)
    positions = []
    for i in range(len(board)):
        for j in range(len(board)):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

positions = find_queen_positions(board)
print(positions)