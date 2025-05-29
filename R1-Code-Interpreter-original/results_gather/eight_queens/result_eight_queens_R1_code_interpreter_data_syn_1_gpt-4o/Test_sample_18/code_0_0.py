def is_valid(board, row, col, n):
    # Check if there's a queen in the same column
    for i in range(row):
        if board[i][col] == 1:
            return False
    
    # Check upper left diagonal
    for i, j in zip(range(row, -1, -1), range(col, -1, -1)):
        if board[i][j] == 1:
            return False
    
    # Check upper right diagonal
    for i, j in zip(range(row, -1, -1), range(col, n)):
        if board[i][j] == 1:
            return False
    
    return True

def solve_n_queens(board, row, n, queens):
    if row == n:
        return True
    
    for col in range(n):
        if board[row][col] == 0 and is_valid(board, row, col, n):
            board[row][col] = 1
            queens.append((row, col))
            if solve_n_queens(board, row + 1, n, queens):
                return True
            board[row][col] = 0
            queens.pop()
    
    return False

def find_queen_positions():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        ['X', 0, 0, 0, 0, 0, 0, 0]
    ]
    
    queens = []
    solve_n_queens(board, 0, n, queens)
    return queens

queens_positions = find_queen_positions()
print(queens_positions)