def is_valid(grid, row, col, num, constraints):
    # Check if num is not in the current row and column
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check constraints
    for (r1, c1, r2, c2, op) in constraints:
        if (r1 == row and c1 == col) or (r2 == row and c2 == col):
            val1 = grid[r1][c1] if grid[r1][c1] != 0 else num
            val2 = grid[r2][c2] if grid[r2][c2] != 0 else num
            if op == '>' and not (val1 > val2):
                return False
            if op == '<' and not (val1 < val2):
                return False
    return True

def solve_futoshiki(grid, constraints, row=0, col=0):
    if row == 6:
        return grid

    if col == 6:
        return solve_futoshiki(grid, constraints, row + 1, 0)

    if grid[row][col] != 0:
        return solve_futoshiki(grid, constraints, row, col + 1)

    for num in range(1, 7):
        if is_valid(grid, row, col, num, constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, constraints, row, col + 1):
                return grid
            grid[row][col] = 0

    return None

# Initial grid with known values
grid = [
    [0, 0, 0, 0, 3, 4],
    [0, 0, 5, 0, 0, 0],
    [5, 2, 0, 0, 0, 3],
    [0, 0, 0, 2, 0, 0],
    [0, 4, 0, 0, 0, 2],
    [0, 0, 0, 0, 2, 0]
]

# Constraints: (row1, col1, row2, col2, operator)
constraints = [
    (0, 2, 0, 3, '>'),
    (1, 1, 1, 2, '>'),
    (2, 1, 2, 2, '>'),
    (2, 4, 2, 5, '>'),
    (3, 0, 3, 1, '<'),
    (3, 3, 3, 4, '<'),
    (3, 4, 3, 5, '<'),
    (4, 5, 4, 4, '>'),
    (5, 0, 5, 1, '>'),
    (5, 1, 5, 2, '<'),
    (5, 4, 5, 5, '<'),
    (0, 0, 1, 0, '<'),
    (1, 2, 2, 2, '<'),
    (2, 5, 3, 5, '<'),
    (3, 1, 4, 1, '<'),
    (4, 5, 5, 5, '<')
]

solution = solve_futoshiki(grid, constraints)
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")