def is_safe(board, row, col):
    # Check this column on upper side
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper diagonal on left side
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False

    # Check upper diagonal on right side
    for i, j in zip(range(row, -1, -1), range(col, len(board))):
        if board[i][j] == 1:
            return False

    return True

def solve_n_queens(board, row):
    if row >= len(board):
        return True

    # If there's already a queen in this row, skip to the next row
    if 1 in board[row]:
        return solve_n_queens(board, row + 1)

    for col in range(len(board)):
        if board[row][col] == 'X':
            continue
        if is_safe(board, row, col):
            board[row][col] = 1
            if solve_n_queens(board, row + 1):
                return True
            board[row][col] = 0

    return False

def find_queen_positions():
    board = [
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 'X', 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    if solve_n_queens(board, 0):
        positions = []
        for i in range(len(board)):
            for j in range(len(board)):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        return positions
    else:
        return "No solution found"

queen_positions = find_queen_positions()
print(queen_positions)