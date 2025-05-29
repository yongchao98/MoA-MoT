def is_safe(board, row, col, n, columns):
    # Check this column on upper side
    if col in columns:
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

def solve_n_queens(board, row, n, columns):
    if row >= n:
        return True

    # If a queen is already placed in this row, move to the next row
    if any(board[row][col] == 1 for col in range(n)):
        return solve_n_queens(board, row + 1, n, columns)

    for col in range(n):
        if board[row][col] == 'X' or board[row][col] == 1:
            continue

        if is_safe(board, row, col, n, columns):
            board[row][col] = 1
            columns.add(col)

            if solve_n_queens(board, row + 1, n, columns):
                return True

            board[row][col] = 0
            columns.remove(col)

    return False

def find_queen_positions(board):
    n = len(board)
    columns = {j for i in range(n) for j in range(n) if board[i][j] == 1}
    solve_n_queens(board, 0, n, columns)
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initial board setup
board = [
    [0, 0, 'X', 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 'X', 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

queen_positions = find_queen_positions(board)
print(queen_positions)