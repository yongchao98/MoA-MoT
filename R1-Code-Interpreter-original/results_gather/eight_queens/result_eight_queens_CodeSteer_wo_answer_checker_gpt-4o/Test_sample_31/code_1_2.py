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

def solve_n_queens(board, row=0, queens_placed=0):
    if queens_placed == 8:
        return True

    if row >= len(board):
        return False

    for col in range(len(board)):
        if board[row][col] == 0 and is_safe(board, row, col):
            board[row][col] = 1
            if solve_n_queens(board, row + 1, queens_placed + 1):
                return True
            board[row][col] = 0

    return False

def find_queen_positions(board):
    # Count already placed queens
    queens_placed = sum(row.count(1) for row in board)
    solve_n_queens(board, queens_placed=queens_placed)
    positions = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initial board setup
board = [
    [0, 0, 0, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# Find and print the positions of the queens
queen_positions = find_queen_positions(board)
print(f"<<<{', '.join(queen_positions)}>>>")