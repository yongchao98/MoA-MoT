def is_valid(board, row, col, queens):
    # Check if the column is free
    for r in range(row):
        if queens[r] == col:
            return False
        # Check diagonals
        if abs(queens[r] - col) == abs(r - row):
            return False
    return True

def solve_n_queens(board, row, queens):
    if row == 8:
        return True
    for col in range(8):
        if board[row][col] == 'X':
            continue
        if is_valid(board, row, col, queens):
            queens[row] = col
            if solve_n_queens(board, row + 1, queens):
                return True
            queens[row] = -1
    return False

def find_queen_positions(board):
    queens = [-1] * 8
    # Place the initial queens
    for r in range(8):
        for c in range(8):
            if board[r][c] == '1':
                queens[r] = c
    solve_n_queens(board, 0, queens)
    return queens

# Initial board setup
board = [
    ['0', '0', 'X', '0', '0', '0', '0', '0'],
    ['1', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', 'X', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0']
]

queen_positions = find_queen_positions(board)
result = ', '.join(f"{r} {c}" for r, c in enumerate(queen_positions))
print(result)