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

    # Skip rows that already have queens
    if any(board[row]):
        return solve_n_queens(board, row + 1, n, positions)

    for i in range(n):
        if is_safe(board, row, i, n):
            board[row][i] = 1
            positions.append((row, i))
            if solve_n_queens(board, row + 1, n, positions):
                return True
            board[row][i] = 0
            positions.pop()

    return False

def find_queen_positions():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    positions = [(1, 0), (3, 5), (6, 2)]
    solve_n_queens(board, 0, n, positions)
    return positions

positions = find_queen_positions()
print(positions)