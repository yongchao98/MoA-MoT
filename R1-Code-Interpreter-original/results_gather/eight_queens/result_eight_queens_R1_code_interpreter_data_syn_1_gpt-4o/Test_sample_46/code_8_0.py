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

def solve_n_queens(board, row, n, queens_positions):
    if row >= n:
        return True

    # Skip rows that already have queens
    if any(board[row][col] == 1 for col in range(n)):
        return solve_n_queens(board, row + 1, n, queens_positions)

    for col in range(n):
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            queens_positions.append((row, col))

            if solve_n_queens(board, row + 1, n, queens_positions):
                return True

            board[row][col] = 0
            queens_positions.pop()

    return False

def find_queen_positions():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, -1],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, 0, 0, 0, 0]
    ]

    # Convert 'X' to a non-placeable marker, e.g., -1
    for i in range(n):
        for j in range(n):
            if board[i][j] == 'X':
                board[i][j] = -1

    queens_positions = [(1, 2), (4, 1)]
    solve_n_queens(board, 0, n, queens_positions)
    return queens_positions

queen_positions = find_queen_positions()
print(queen_positions)