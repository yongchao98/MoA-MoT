def is_safe(board, row, col):
    # Check if the position is marked with 'X'
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

            # Backtrack
            board[row][col] = 0
            positions.pop()

    return False

def place_queens(initial_board):
    board = [list(row) for row in initial_board]
    positions = []

    # Add pre-placed queens to positions
    for r in range(len(board)):
        for c in range(len(board[r])):
            if board[r][c] == 1:
                positions.append((r, c))

    # Find the first empty row to start placing queens
    start_row = next((r for r in range(len(board)) if 1 not in board[r]), 0)

    if solve_n_queens(board, start_row, positions):
        result = ', '.join(f'{r} {c}' for r, c in positions)
        print(f'<<<{result}>>>')
    else:
        print("No solution exists")

# Example initial board setup
initial_board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

place_queens(initial_board)