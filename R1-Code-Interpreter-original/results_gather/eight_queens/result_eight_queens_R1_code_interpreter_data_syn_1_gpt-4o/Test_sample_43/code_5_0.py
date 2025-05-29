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

def solve_n_queens(board, row, n, result):
    if row >= n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(board, row, col, n):
            board[row][col] = 1
            result.append((row, col))
            if solve_n_queens(board, row + 1, n, result):
                return True
            board[row][col] = 0
            result.pop()

    return False

def find_queen_positions():
    board = [
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 'X'],
        [0, 0, 'X', 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0]
    ]

    n = 8
    result = []

    # Mark the existing queen positions
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                result.append((i, j))

    # Start solving from the first row that doesn't have a queen
    start_row = 0
    for i in range(n):
        if any(board[i][j] == 1 for j in range(n)):
            start_row = i + 1

    # Solve the problem starting from the first empty row
    solve_n_queens(board, start_row, n, result)

    return result

queen_positions = find_queen_positions()
print(queen_positions)