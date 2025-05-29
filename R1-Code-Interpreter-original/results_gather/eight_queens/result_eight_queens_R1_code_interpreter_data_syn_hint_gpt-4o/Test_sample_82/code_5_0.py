def is_safe(row, col, cols, diag1, diag2):
    return col not in cols and (row - col) not in diag1 and (row + col) not in diag2

def solve_n_queens(board, row, n, cols, diag1, diag2):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)

            if solve_n_queens(board, row + 1, n, cols, diag1, diag2):
                return True

            board[row][col] = 0
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)

    return False

def find_queen_positions(board):
    n = len(board)
    cols = set()
    diag1 = set()
    diag2 = set()

    # Mark the initial queen's position
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                cols.add(j)
                diag1.add(i - j)
                diag2.add(i + j)

    solve_n_queens(board, 0, n, cols, diag1, diag2)
    positions = []
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    return positions

# Initial board setup
board = [
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, -1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, -1, 0]
]

# Convert 'X' to a non-placeable value
for i in range(len(board)):
    for j in range(len(board[i])):
        if board[i][j] == -1:
            board[i][j] = -1

# Find and print the positions of the queens
queen_positions = find_queen_positions(board)
print(queen_positions)