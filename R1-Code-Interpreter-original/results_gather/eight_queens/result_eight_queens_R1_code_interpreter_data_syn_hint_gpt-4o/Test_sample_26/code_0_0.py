def is_valid(board, row, col, n):
    # Check if there's a queen in the same column
    for i in range(row):
        if board[i][col] == 1:
            return False

    # Check upper left diagonal
    i, j = row, col
    while i >= 0 and j >= 0:
        if board[i][j] == 1:
            return False
        i -= 1
        j -= 1

    # Check upper right diagonal
    i, j = row, col
    while i >= 0 and j < n:
        if board[i][j] == 1:
            return False
        i -= 1
        j += 1

    return True

def solve_n_queens(board, row, n, positions):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_valid(board, row, col, n):
            board[row][col] = 1
            positions.append((row, col))
            if solve_n_queens(board, row + 1, n, positions):
                return True
            board[row][col] = 0
            positions.pop()

    return False

def find_queen_positions():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0]
    ]

    positions = [(1, 0), (7, 1)]
    solve_n_queens(board, 0, n, positions)
    return positions

queen_positions = find_queen_positions()
print(queen_positions)