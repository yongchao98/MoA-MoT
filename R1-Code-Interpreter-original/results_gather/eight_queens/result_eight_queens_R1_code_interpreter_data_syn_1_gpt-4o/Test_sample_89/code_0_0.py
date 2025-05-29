def is_safe(board, row, col, queens):
    # Check if the column is safe
    for r, c in queens:
        if c == col or abs(r - row) == abs(c - col):
            return False
    return True

def solve_n_queens(board, row, queens):
    if row == 8:
        return queens
    
    for col in range(8):
        if board[row][col] == 'X':
            continue
        if is_safe(board, row, col, queens):
            queens.append((row, col))
            result = solve_n_queens(board, row + 1, queens)
            if result:
                return result
            queens.pop()
    return None

# Initial board setup
board = [
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', 'X', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', 'X'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', 'X'],
    ['0', 'X', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0']
]

# Solve the problem
queens = solve_n_queens(board, 0, [])
print(queens)