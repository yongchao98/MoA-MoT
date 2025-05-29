def is_safe(row, col, columns, diag1, diag2, blocked):
    if col in columns or (row - col) in diag1 or (row + col) in diag2 or (row, col) in blocked:
        return False
    return True

def solve_queens(n, row, columns, diag1, diag2, blocked, queens_positions):
    if row == n:
        return True

    for col in range(n):
        if is_safe(row, col, columns, diag1, diag2, blocked):
            columns.add(col)
            diag1.add(row - col)
            diag2.add(row + col)
            queens_positions.append((row, col))

            if solve_queens(n, row + 1, columns, diag1, diag2, blocked, queens_positions):
                return True

            columns.remove(col)
            diag1.remove(row - col)
            diag2.remove(row + col)
            queens_positions.pop()

    return False

def find_queen_positions():
    n = 8
    blocked = {(1, 1), (2, 4), (7, 2)}
    queens_positions = [(3, 6)]
    columns = {6}
    diag1 = {-3}
    diag2 = {9}

    # Start solving from row 0, but skip row 3 since it already has a queen
    for row in range(n):
        if row == 3:
            continue
        if solve_queens(n, row, columns, diag1, diag2, blocked, queens_positions):
            break

    return queens_positions

queen_positions = find_queen_positions()
print(queen_positions)