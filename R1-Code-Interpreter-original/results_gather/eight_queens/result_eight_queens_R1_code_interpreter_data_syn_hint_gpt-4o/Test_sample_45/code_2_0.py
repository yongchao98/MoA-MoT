def is_safe(row, col, columns, diag1, diag2, board):
    if col in columns or (row - col) in diag1 or (row + col) in diag2 or board[row][col] == 'X':
        return False
    return True

def solve_n_queens(board, row, n, columns, diag1, diag2, positions):
    if row == n:
        return True

    for col in range(n):
        if is_safe(row, col, columns, diag1, diag2, board):
            board[row][col] = 1
            columns.add(col)
            diag1.add(row - col)
            diag2.add(row + col)
            positions.append((row, col))

            if solve_n_queens(board, row + 1, n, columns, diag1, diag2, positions):
                return True

            # Backtrack
            board[row][col] = 0
            columns.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)
            positions.pop()

    return False

def find_queen_positions(board):
    n = len(board)
    columns = set()
    diag1 = set()
    diag2 = set()
    positions = []

    # Pre-fill the board with existing queens
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                columns.add(j)
                diag1.add(i - j)
                diag2.add(i + j)
                positions.append((i, j))

    solve_n_queens(board, 0, n, columns, diag1, diag2, positions)
    return [f"{i} {j}" for i, j in positions]

# Initial board setup
board = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 'X', 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 'X', 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 'X', 0, 'X', 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

queen_positions = find_queen_positions(board)
print(queen_positions)