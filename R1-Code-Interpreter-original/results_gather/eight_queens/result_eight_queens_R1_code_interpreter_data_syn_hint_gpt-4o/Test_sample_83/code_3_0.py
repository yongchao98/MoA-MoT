def is_safe(row, col, cols, diag1, diag2):
    return not (cols[col] or diag1[row - col] or diag2[row + col])

def place_queens(board, row, cols, diag1, diag2, positions):
    n = len(board)
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            cols[col] = diag1[row - col] = diag2[row + col] = True
            positions.append((row, col))

            if place_queens(board, row + 1, cols, diag1, diag2, positions):
                return True

            board[row][col] = 0
            cols[col] = diag1[row - col] = diag2[row + col] = False
            positions.pop()

    return False

def find_queen_positions():
    n = 8
    board = [[0] * n for _ in range(n)]
    board[0][4] = 1
    board[4][7] = 1

    cols = [False] * n
    diag1 = [False] * (2 * n - 1)
    diag2 = [False] * (2 * n - 1)

    cols[4] = True
    cols[7] = True
    diag1[0 - 4] = True
    diag1[4 - 7] = True
    diag2[0 + 4] = True
    diag2[4 + 7] = True

    positions = [(0, 4), (4, 7)]

    place_queens(board, 1, cols, diag1, diag2, positions)

    return [f"{r} {c}" for r, c in positions]

print(find_queen_positions())