def is_safe(row, col, cols, diag1, diag2):
    return not (col in cols or (row - col) in diag1 or (row + col) in diag2)

def place_queens(n, row, cols, diag1, diag2, positions):
    if row == n:
        return True

    for col in range(n):
        if is_safe(row, col, cols, diag1, diag2):
            positions.append((row, col))
            cols.add(col)
            diag1.add(row - col)
            diag2.add(row + col)

            if place_queens(n, row + 1, cols, diag1, diag2, positions):
                return True

            positions.pop()
            cols.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)

    return False

def find_queen_positions():
    n = 8
    positions = [(0, 4), (4, 7)]
    cols = {4, 7}
    diag1 = {0 - 4, 4 - 7}
    diag2 = {0 + 4, 4 + 7}

    place_queens(n, 1, cols, diag1, diag2, positions)

    return [f"{r} {c}" for r, c in positions]

print(find_queen_positions())