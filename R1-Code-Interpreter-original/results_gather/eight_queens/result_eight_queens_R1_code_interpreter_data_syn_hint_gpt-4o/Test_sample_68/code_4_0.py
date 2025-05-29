def is_safe(queens, row, col):
    for r in range(row):
        c = queens[r]
        if c == col or abs(c - col) == abs(r - row):
            return False
    return True

def solve_queens(queens, row, n):
    if row == n:
        return True

    if row == 4:  # Skip the row with the initial queen
        return solve_queens(queens, row + 1, n)

    for col in range(n):
        if (row, col) == (5, 5):  # Skip restricted position
            continue
        if is_safe(queens, row, col):
            queens[row] = col
            if solve_queens(queens, row + 1, n):
                return True
            queens[row] = -1

    return False

def find_queen_positions():
    n = 8
    queens = [-1] * n
    queens[4] = 2  # Initial placement of the queen

    solve_queens(queens, 0, n)

    positions = []
    for i in range(n):
        if queens[i] != -1:
            positions.append(f"{i} {queens[i]}")

    return positions

print(find_queen_positions())