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

def solve_n_queens(board, row, n, queens_count):
    if queens_count == 8:
        return True

    for col in range(n):
        if board[row][col] == 'X' or board[row][col] == 1:
            continue
        if is_safe(board, row, col, n):
            board[row][col] = 1
            if solve_n_queens(board, row + 1, n, queens_count + 1):
                return True
            board[row][col] = 0

    return False

def find_queen_positions(board):
    n = len(board)
    # Count the initial number of queens
    initial_queens = sum(row.count(1) for row in board)
    # Start solving from the first row
    solve_n_queens(board, 0, n, initial_queens)
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    ['X', 0, 'X', 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# Find and print the positions of the queens
queen_positions = find_queen_positions(board)
print(f"<<<{', '.join(queen_positions)}>>>")