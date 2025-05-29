def is_safe(row, col, cols, diag1, diag2):
    return not (cols[col] or diag1[row - col] or diag2[row + col])

def place_queens(board, row, cols, diag1, diag2, positions, pre_placed_rows):
    n = len(board)
    if row == n:
        return True

    # Skip rows that already have pre-placed queens
    if row in pre_placed_rows:
        return place_queens(board, row + 1, cols, diag1, diag2, positions, pre_placed_rows)

    for col in range(n):
        if board[row][col] == 0 and is_safe(row, col, cols, diag1, diag2):
            board[row][col] = 1
            cols[col] = diag1[row - col] = diag2[row + col] = True
            positions.append(f"{row} {col}")

            if place_queens(board, row + 1, cols, diag1, diag2, positions, pre_placed_rows):
                return True

            board[row][col] = 0
            cols[col] = diag1[row - col] = diag2[row + col] = False
            positions.pop()

    return False

def find_queen_positions():
    n = 8
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

    cols = [False] * n
    diag1 = [False] * (2 * n)
    diag2 = [False] * (2 * n)

    # Mark the pre-placed queens
    pre_placed = [(1, 2), (5, 4), (6, 0)]
    pre_placed_rows = {r for r, c in pre_placed}
    for r, c in pre_placed:
        board[r][c] = 1
        cols[c] = diag1[r - c] = diag2[r + c] = True

    positions = [f"{r} {c}" for r, c in pre_placed]

    # Start placing queens from the first row
    place_queens(board, 0, cols, diag1, diag2, positions, pre_placed_rows)

    return positions

print(find_queen_positions())