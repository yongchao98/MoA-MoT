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

def solve_n_queens(board, row, positions):
    if row >= len(board):
        return True

    # Get columns sorted by the number of available positions
    available_cols = [col for col in range(len(board)) if board[row][col] == 0]
    available_cols.sort(key=lambda col: sum(board[i][col] == 0 for i in range(row, len(board))))

    for col in available_cols:
        if is_safe(board, row, col):
            board[row][col] = 1
            positions.append((row, col))
            if solve_n_queens(board, row + 1, positions):
                return True
            board[row][col] = 0
            positions.pop()

    return False

def place_queens():
    # Initialize the board with 0s and -1s for blocked positions
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    # Example of blocked positions (X)
    blocked_positions = [(0, 1), (1, 3), (2, 5)]  # Add actual blocked positions here
    for r, c in blocked_positions:
        board[r][c] = -1

    positions = []

    if solve_n_queens(board, 0, positions):
        result = ', '.join(f'{r} {c}' for r, c in positions)
        print(f'<<<{result}>>>')
    else:
        print("No solution exists")

place_queens()