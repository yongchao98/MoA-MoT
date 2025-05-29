def is_safe(board, row, col):
    # Check this row on left side
    for i in range(col):
        if board[row][i] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check lower diagonal on left side
    for i, j in zip(range(row, len(board), 1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, col, positions):
    if col >= len(board):
        return True

    # If there's already a queen in this column, move to the next column
    for i in range(len(board)):
        if board[i][col] == 1:
            return solve_n_queens(board, col + 1, positions)

    for i in range(len(board)):
        if board[i][col] == 0 and is_safe(board, i, col):
            board[i][col] = 1
            positions.append((i, col))

            if solve_n_queens(board, col + 1, positions):
                return True

            board[i][col] = 0
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

    positions = [(0, 2)]
    # Start from the first column that doesn't have a queen
    solve_n_queens(board, 1, positions)
    return positions

queen_positions = find_queen_positions()
print(queen_positions)