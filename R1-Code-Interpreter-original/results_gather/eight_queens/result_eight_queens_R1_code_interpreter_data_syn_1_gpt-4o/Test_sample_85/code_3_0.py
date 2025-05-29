def is_safe(board, row, col, rows, cols, diag1, diag2):
    if rows[row] or cols[col] or diag1[row - col] or diag2[row + col]:
        return False
    return True

def solve_n_queens(board, col, positions, rows, cols, diag1, diag2):
    if col >= len(board):
        return True

    for row in range(len(board)):
        if board[row][col] == 0 and is_safe(board, row, col, rows, cols, diag1, diag2):
            board[row][col] = 1
            positions.append((row, col))
            rows[row] = cols[col] = diag1[row - col] = diag2[row + col] = True

            if solve_n_queens(board, col + 1, positions, rows, cols, diag1, diag2):
                return True

            board[row][col] = 0
            positions.pop()
            rows[row] = cols[col] = diag1[row - col] = diag2[row + col] = False

    return False

def find_queen_positions():
    board = [
        [0, 'X', 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 'X', 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 'X', 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    positions = [(0, 2)]
    rows = [False] * 8
    cols = [False] * 8
    diag1 = [False] * 15
    diag2 = [False] * 15

    # Mark the initial queen positions
    for r, c in positions:
        rows[r] = cols[c] = diag1[r - c] = diag2[r + c] = True

    solve_n_queens(board, 0, positions, rows, cols, diag1, diag2)
    return positions

queen_positions = find_queen_positions()
print(queen_positions)