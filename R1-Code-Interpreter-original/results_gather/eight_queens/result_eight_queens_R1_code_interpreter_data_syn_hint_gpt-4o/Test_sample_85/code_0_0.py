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

def solve_n_queens(board, row, n, positions):
    if row >= n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            positions.append((row, col))

            if solve_n_queens(board, row + 1, n, positions):
                return True

            # Backtrack
            board[row][col] = 0
            positions.pop()

    return False

def find_queen_positions():
    board = [
        [0, 'X', 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 'X', 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    n = 8
    positions = []

    # Mark initial queens
    initial_queens = [(0, 2)]
    for r, c in initial_queens:
        board[r][c] = 1
        positions.append((r, c))

    if solve_n_queens(board, 0, n, positions):
        return positions
    else:
        return "No solution found"

positions = find_queen_positions()
print(positions)