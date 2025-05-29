def is_safe(board, row, col, occupied_columns):
    # Check if the column is already occupied
    if col in occupied_columns:
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

def solve_n_queens(board, row, positions, occupied_columns):
    if row >= len(board):
        return True

    # Skip the row with the pre-placed queen
    if row == 4:
        return solve_n_queens(board, row + 1, positions, occupied_columns)

    for col in range(len(board)):
        if board[row][col] == 0 and is_safe(board, row, col, occupied_columns):
            board[row][col] = 1
            positions.append((row, col))
            occupied_columns.add(col)
            if solve_n_queens(board, row + 1, positions, occupied_columns):
                return True
            board[row][col] = 0
            positions.pop()
            occupied_columns.remove(col)

    return False

def place_queens():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    positions = [(4, 5)]
    occupied_columns = {5}  # Column 5 is occupied by the pre-placed queen

    if solve_n_queens(board, 0, positions, occupied_columns):
        print("<<<" + ", ".join(f"{r} {c}" for r, c in positions) + ">>>")
    else:
        print("No solution found")

place_queens()