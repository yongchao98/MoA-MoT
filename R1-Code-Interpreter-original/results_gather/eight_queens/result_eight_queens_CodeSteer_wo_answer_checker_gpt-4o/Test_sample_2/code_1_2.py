def is_safe(board, row, col):
    # Check if the position is restricted
    if board[row][col] == 'X':
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
            board[row][col] = 0
            positions.pop()

    return False

def place_queens():
    # Initialize the board with given conditions
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

    # Example of pre-placed queens and restricted positions
    # board[0][0] = 1  # Pre-placed queen
    # board[1][2] = 'X'  # Restricted position

    positions = [(i, j) for i in range(8) for j in range(8) if board[i][j] == 1]

    if solve_n_queens(board, 0, positions):
        result = ', '.join(f'{r} {c}' for r, c in positions)
        print(f'<<<{result}>>>')
    else:
        print("No solution exists")

place_queens()