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

def solve_n_queens(board, row, queens_placed):
    if queens_placed == 8:
        return True

    if row >= len(board):
        return False

    for col in range(len(board)):
        if board[row][col] == 'X' or board[row][col] == 1:
            continue
        if is_safe(board, row, col):
            board[row][col] = 1
            if solve_n_queens(board, row + 1, queens_placed + 1):
                return True
            board[row][col] = 0

    return solve_n_queens(board, row + 1, queens_placed)

def find_queen_positions():
    board = [[0 for _ in range(8)] for _ in range(8)]
    board[4][0] = 1  # Pre-placed queen
    board[7][5] = 'X'  # Blocked position

    if solve_n_queens(board, 0, 1):  # Start with one queen already placed
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