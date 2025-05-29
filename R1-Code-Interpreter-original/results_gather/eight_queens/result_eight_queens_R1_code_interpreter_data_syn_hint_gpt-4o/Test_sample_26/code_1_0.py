def is_valid(board, row, col, n, cols, diag1, diag2):
    if cols[col] or diag1[row - col] or diag2[row + col]:
        return False
    return True

def solve_n_queens(board, row, n, positions, cols, diag1, diag2):
    if len(positions) == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_valid(board, row, col, n, cols, diag1, diag2):
            board[row][col] = 1
            positions.append((row, col))
            cols[col] = diag1[row - col] = diag2[row + col] = True

            if solve_n_queens(board, row + 1, n, positions, cols, diag1, diag2):
                return True

            board[row][col] = 0
            positions.pop()
            cols[col] = diag1[row - col] = diag2[row + col] = False

    return False

def find_queen_positions():
    n = 8
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0]
    ]

    positions = [(1, 0), (7, 1)]
    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)

    # Mark initial queens
    for r, c in positions:
        cols[c] = diag1[r - c] = diag2[r + c] = True

    # Start solving from the first row
    solve_n_queens(board, 0, n, positions, cols, diag1, diag2)
    return positions

queen_positions = find_queen_positions()
print(queen_positions)