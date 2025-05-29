def is_safe(row, col, cols, diag1, diag2):
    return not (cols[col] or diag1[row - col] or diag2[row + col])

def place_queens(n, row, board, cols, diag1, diag2, positions):
    if row == n:
        return True

    if row == 2:  # Skip row 2 as it already has a queen
        return place_queens(n, row + 1, board, cols, diag1, diag2, positions)

    for col in range(n):
        if is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            cols[col] = diag1[row - col] = diag2[row + col] = True
            positions.append((row, col))

            if place_queens(n, row + 1, board, cols, diag1, diag2, positions):
                return True

            board[row][col] = 0
            cols[col] = diag1[row - col] = diag2[row + col] = False
            positions.pop()

    return False

# Initialize the board and constraints
n = 8
board = [[0] * n for _ in range(n)]
cols = [False] * n
diag1 = [False] * (2 * n)
diag2 = [False] * (2 * n)
positions = []

# Place the initial queen
board[2][3] = 1
cols[3] = diag1[2 - 3] = diag2[2 + 3] = True
positions.append((2, 3))

# Solve the problem
place_queens(n, 0, board, cols, diag1, diag2, positions)

print(positions)