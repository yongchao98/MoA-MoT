def is_safe(board, row, col):
    # Check this column on upper side
    for i in range(len(board)):
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

    # Skip the row where the queen is already placed
    if row == 2:
        return solve_n_queens(board, row + 1)

    for col in range(len(board)):
        if is_safe(board, row, col):
            board[row][col] = 1
            if solve_n_queens(board, row + 1):
                return True
            board[row][col] = 0

    return False

def find_queen_positions():
    board = [[0] * 8 for _ in range(8)]
    board[2][6] = 1  # Pre-placed queen

    if solve_n_queens(board, 0):
        positions = []
        for i in range(8):
            for j in range(8):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        return positions
    else:
        return "No solution found"

queen_positions = find_queen_positions()
print(queen_positions)