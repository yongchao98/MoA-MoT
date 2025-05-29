def is_safe(board, row, col):
    # Check if the cell is restricted
    if board[row][col] == -1:
        return False

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

def solve_n_queens(board, row, positions):
    if row >= len(board):
        return True

    for col in range(len(board)):
        if board[row][col] == 0 and is_safe(board, row, col):
            board[row][col] = 1
            positions.append((row, col))

            if solve_n_queens(board, row + 1, positions):
                return True

            # Backtrack
            board[row][col] = 0
            positions.pop()

    return False

def place_queens():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0],  # 'X' is represented by -1
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    positions = []
    if solve_n_queens(board, 0, positions):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
    else:
        print("No solution found")

place_queens()