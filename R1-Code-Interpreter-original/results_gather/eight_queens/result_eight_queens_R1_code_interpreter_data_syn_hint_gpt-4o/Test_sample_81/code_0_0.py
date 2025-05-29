def is_valid(board, row, col):
    # Check column
    for i in range(row):
        if board[i][col] == 1:
            return False
    # Check diagonal (top-left to bottom-right)
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    # Check diagonal (top-right to bottom-left)
    for i, j in zip(range(row, -1, -1), range(col, len(board))):
        if board[i][j] == 1:
            return False
    return True

def solve_n_queens(board, row, queens_positions):
    if row == len(board):
        return True
    for col in range(len(board)):
        if board[row][col] == 0 and is_valid(board, row, col):
            board[row][col] = 1
            queens_positions.append((row, col))
            if solve_n_queens(board, row + 1, queens_positions):
                return True
            board[row][col] = 0
            queens_positions.pop()
    return False

def find_queen_positions():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 'X', 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 'X', 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]
    queens_positions = []
    solve_n_queens(board, 0, queens_positions)
    return queens_positions

queen_positions = find_queen_positions()
print(queen_positions)