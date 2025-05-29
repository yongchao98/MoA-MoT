def is_valid(board, row, col):
    # Check column
    for i in range(row):
        if board[i] == col:
            return False
    # Check diagonals
    for i in range(row):
        if abs(board[i] - col) == abs(i - row):
            return False
    return True

def solve_n_queens(board, row, positions):
    if row == 8:
        return True
    for col in range(8):
        if board[row][col] == 'X':
            continue
        if is_valid(positions, row, col):
            positions[row] = col
            if solve_n_queens(board, row + 1, positions):
                return True
            positions[row] = -1
    return False

def place_queens():
    board = [
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]
    positions = [-1] * 8
    positions[0] = 1  # Queen already placed at (0, 1)
    if solve_n_queens(board, 1, positions):
        result = ', '.join(f'{i} {positions[i]}' for i in range(8))
        print(f'<<<{result}>>>')
    else:
        print("No solution found")

place_queens()