def is_safe(row, col, cols, diag1, diag2):
    return not (cols[col] or diag1[row - col] or diag2[row + col])

def place_queens(board, row, n, cols, diag1, diag2, positions):
    if row == n:
        return True

    for col in range(n):
        if board[row][col] == 0 and is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            cols[col] = diag1[row - col] = diag2[row + col] = True
            positions.append((row, col))

            if place_queens(board, row + 1, n, cols, diag1, diag2, positions):
                return True

            # Backtrack
            board[row][col] = 0
            cols[col] = diag1[row - col] = diag2[row + col] = False
            positions.pop()

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

    n = 8
    positions = []

    # Initialize tracking arrays
    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)

    # Mark initial queens
    initial_queens = [(0, 2)]
    for r, c in initial_queens:
        board[r][c] = 1
        positions.append((r, c))
        cols[c] = diag1[r - c] = diag2[r + c] = True

    if place_queens(board, 0, n, cols, diag1, diag2, positions):
        return positions
    else:
        return "No solution found"

positions = find_queen_positions()
print(positions)