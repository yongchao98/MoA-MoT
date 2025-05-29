def is_safe(board, row, col, queens):
    # Check this column on upper side
    for r, c in queens:
        if c == col or abs(r - row) == abs(c - col):
            return False
    return True

def solve_n_queens(board, row, queens):
    if len(queens) == 8:
        return True

    if row >= len(board):
        return False

    for col in range(len(board)):
        if board[row][col] == 0 and is_safe(board, row, col, queens):
            queens.append((row, col))
            if solve_n_queens(board, row + 1, queens):
                return True
            queens.pop()

    return solve_n_queens(board, row + 1, queens)

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 'X', 0, 0],
    [0, 0, 'X', 0, 0, 0, 0, 0]
]

# Start with the queens already placed
queens = [(4, 4), (6, 0)]

# Solve the problem
solve_n_queens(board, 0, queens)

# Collect the positions of the queens
queen_positions = [f"{r} {c}" for r, c in queens]

print(queen_positions)