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

def solve_n_queens(board, col, n):
    if col >= n:
        return True

    for i in range(n):
        if board[i][col] == 0 and is_safe(board, i, col, n):
            board[i][col] = 1
            if solve_n_queens(board, col + 1, n):
                return True
            board[i][col] = 0

    return False

def find_queen_positions(board):
    n = len(board)
    # Mark pre-placed queens
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                if not is_safe(board, i, j, n):
                    return "No solution"
    
    # Start solving from the first column
    if solve_n_queens(board, 0, n):
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        return positions
    else:
        return "No solution"

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, -1],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [-1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, -1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# Find the positions of the queens
queen_positions = find_queen_positions(board)
print(queen_positions)