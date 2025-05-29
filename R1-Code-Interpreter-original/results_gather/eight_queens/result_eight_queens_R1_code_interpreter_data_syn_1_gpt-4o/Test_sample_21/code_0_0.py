def is_safe(board, row, col, queens):
    # Check if the column is safe
    for i in range(row):
        if queens[i] == col:
            return False
        # Check diagonals
        if abs(queens[i] - col) == abs(i - row):
            return False
    return True

def solve_n_queens(board, row, queens):
    if row == 8:
        return True
    
    for col in range(8):
        if board[row][col] == 'X':
            continue
        if is_safe(board, row, col, queens):
            queens[row] = col
            if solve_n_queens(board, row + 1, queens):
                return True
            queens[row] = -1
    return False

# Initial board setup
board = [
    ['0', '0', '0', '1', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', 'X', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0']
]

# Initial queen positions
queens = [-1] * 8
queens[0] = 3  # Queen already placed at (0, 3)

# Solve the problem
solve_n_queens(board, 1, queens)

# Prepare the result
result = []
for i in range(8):
    result.append(f"{i} {queens[i]}")

print(f"<<<{', '.join(result)}>>>")