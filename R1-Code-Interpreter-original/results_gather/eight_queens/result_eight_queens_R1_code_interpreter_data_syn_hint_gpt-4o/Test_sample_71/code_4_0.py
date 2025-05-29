def is_safe(row, col, cols, diag1, diag2):
    return not (cols[col] or diag1[row - col] or diag2[row + col])

def solve_n_queens(n, board, row, cols, diag1, diag2, positions):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == -1:
            continue
        if is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            cols[col] = diag1[row - col] = diag2[row + col] = True
            positions.append((row, col))
            if solve_n_queens(n, board, row + 1, cols, diag1, diag2, positions):
                return True
            board[row][col] = 0
            cols[col] = diag1[row - col] = diag2[row + col] = False
            positions.pop()

    return False

def main():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -1, 0]
    ]

    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)
    positions = []

    solve_n_queens(n, board, 0, cols, diag1, diag2, positions)
    print(positions)

main()