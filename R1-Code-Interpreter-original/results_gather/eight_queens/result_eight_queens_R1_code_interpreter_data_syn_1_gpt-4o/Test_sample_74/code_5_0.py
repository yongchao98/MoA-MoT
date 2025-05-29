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

def solve_n_queens(board, row, n):
    if row >= n:
        return True

    # If a queen is already placed in this row, move to the next row
    if any(board[row][i] == 1 for i in range(n)):
        return solve_n_queens(board, row + 1, n)

    for i in range(n):
        if board[row][i] == 0 and is_safe(board, row, i, n):
            board[row][i] = 1
            if solve_n_queens(board, row + 1, n):
                return True
            board[row][i] = 0

    return False

def find_queen_positions(board):
    n = len(board)
    solve_n_queens(board, 0, n)
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append((i, j))
    return positions

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0]
]

positions = find_queen_positions(board)
print(positions)