def is_safe(row, col, cols, diag1, diag2, board):
    if col in cols or (row - col) in diag1 or (row + col) in diag2 or board[row][col] == 'X':
        return False
    return True

def solve_n_queens(board, row, n, cols, diag1, diag2, positions):
    if row == n:
        return True

    for col in range(n):
        if is_safe(row, col, cols, diag1, diag2, board):
            board[row][col] = 1
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)
            positions.append((row, col))

            if solve_n_queens(board, row + 1, n, cols, diag1, diag2, positions):
                return True

            # Backtrack
            board[row][col] = 0
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)
            positions.pop()

    return False

def find_queen_positions(board):
    n = len(board)
    cols = set()
    diag1 = set()
    diag2 = set()
    positions = []

    # Mark existing queens
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                cols.add(j)
                diag1.add(i - j)
                diag2.add(i + j)
                positions.append((i, j))

    # Start solving from the first row
    solve_n_queens(board, 0, n, cols, diag1, diag2, positions)
    return positions

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

positions = find_queen_positions(board)
formatted_positions = ', '.join(f"{r} {c}" for r, c in positions)
print(f"<<<{formatted_positions}>>>")