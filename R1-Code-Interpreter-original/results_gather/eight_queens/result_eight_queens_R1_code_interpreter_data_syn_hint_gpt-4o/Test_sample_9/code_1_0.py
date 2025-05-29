def is_valid(board, row, col):
    # Check if the position is marked with 'X'
    if board[row][col] == 'X':
        return False
    # Check column
    for i in range(row):
        if board[i][col] == 1:
            return False
    # Check left diagonal
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    # Check right diagonal
    for i, j in zip(range(row, -1, -1), range(col, len(board))):
        if board[i][j] == 1:
            return False
    return True

def solve_n_queens(board, row, positions):
    if row == len(board):
        return True
    for col in range(len(board)):
        if is_valid(board, row, col):
            board[row][col] = 1
            positions.append((row, col))
            if solve_n_queens(board, row + 1, positions):
                return True
            board[row][col] = 0
            positions.pop()
    return False

def find_queen_positions():
    # Initialize the board with '0' and 'X' for blocked positions
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
    positions = []
    solve_n_queens(board, 0, positions)
    return positions

queen_positions = find_queen_positions()
print(queen_positions)