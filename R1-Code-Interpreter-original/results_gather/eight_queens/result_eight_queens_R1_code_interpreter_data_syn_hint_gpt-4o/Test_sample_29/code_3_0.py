def is_safe(row, col, columns, main_diag, anti_diag):
    if columns[col] or main_diag[row - col] or anti_diag[row + col]:
        return False
    return True

def place_queens(n, row, board, columns, main_diag, anti_diag, positions):
    if row == n:
        return True

    for col in range(n):
        if is_safe(row, col, columns, main_diag, anti_diag):
            board[row][col] = 1
            columns[col] = True
            main_diag[row - col] = True
            anti_diag[row + col] = True
            positions.append((row, col))

            if place_queens(n, row + 1, board, columns, main_diag, anti_diag, positions):
                return True

            board[row][col] = 0
            columns[col] = False
            main_diag[row - col] = False
            anti_diag[row + col] = False
            positions.pop()

    return False

def find_queen_positions():
    n = 8
    board = [[0] * n for _ in range(n)]
    positions = [(1, 0), (3, 5), (6, 2)]
    columns = [False] * n
    main_diag = [False] * (2 * n)
    anti_diag = [False] * (2 * n)

    # Mark initial positions
    for r, c in positions:
        board[r][c] = 1
        columns[c] = True
        main_diag[r - c] = True
        anti_diag[r + c] = True

    # Start placing queens from the first unoccupied row
    place_queens(n, 0, board, columns, main_diag, anti_diag, positions)
    return positions

positions = find_queen_positions()
print(positions)