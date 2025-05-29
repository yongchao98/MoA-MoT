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
        if board[row][col] == 'X':
            continue
        if is_safe(board, row, col, n):
            board[row][col] = 1
            positions.append((row, col))
            if solve_n_queens(board, row + 1, n, positions):
                return True
            board[row][col] = 0
            positions.pop()

    return False

def find_queen_positions():
    n = 8
    board = [[0 for _ in range(n)] for _ in range(n)]
    board[5][4] = 'X'  # Mark the cell with 'X'
    positions = []
    solve_n_queens(board, 0, n, positions)
    return positions

queen_positions = find_queen_positions()
print(queen_positions)