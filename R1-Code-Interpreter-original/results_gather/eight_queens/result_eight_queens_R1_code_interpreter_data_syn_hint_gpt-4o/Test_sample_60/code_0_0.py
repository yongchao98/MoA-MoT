def is_valid(board, row, col):
    # Check column
    for i in range(row):
        if board[i] == col:
            return False
    # Check diagonal
    for i in range(row):
        if abs(board[i] - col) == abs(i - row):
            return False
    return True

def solve_n_queens(board, row):
    if row == len(board):
        return True
    for col in range(len(board)):
        if is_valid(board, row, col):
            board[row] = col
            if solve_n_queens(board, row + 1):
                return True
            board[row] = -1
    return False

def find_queen_positions():
    board = [-1] * 8
    board[7] = 7  # Queen already placed at (7,7)
    solve_n_queens(board, 0)
    return [(i, board[i]) for i in range(8)]

queen_positions = find_queen_positions()
print(queen_positions)