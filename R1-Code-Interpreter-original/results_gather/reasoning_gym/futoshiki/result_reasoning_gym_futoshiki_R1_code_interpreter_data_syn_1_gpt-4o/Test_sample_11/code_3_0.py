def is_valid(grid, row, col, num, constraints):
    # Check if num is not in the current row and column
    for i in range(4):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check horizontal constraints
    for (r, c1, c2, op) in constraints['horizontal']:
        if r == row:
            if c1 == col:
                if op == '<' and not (num < grid[r][c2]):
                    return False
                if op == '>' and not (num > grid[r][c2]):
                    return False
            if c2 == col:
                if op == '<' and not (grid[r][c1] < num):
                    return False
                if op == '>' and not (grid[r][c1] > num):
                    return False

    # Check vertical constraints
    for (c, r1, r2, op) in constraints['vertical']:
        if c == col:
            if r1 == row:
                if op == '∧' and not (num < grid[r2][c]):
                    return False
                if op == '∨' and not (num > grid[r2][c]):
                    return False
            if r2 == row:
                if op == '∧' and not (grid[r1][c] < num):
                    return False
                if op == '∨' and not (grid[r1][c] > num):
                    return False

    return True

def solve_futoshiki(grid, row, col, constraints):
    if row == 4:
        return True
    if col == 4:
        return solve_futoshiki(grid, row + 1, 0, constraints)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1, constraints)

    for num in range(1, 5):
        if is_valid(grid, row, col, num, constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1, constraints):
                return True
            grid[row][col] = 0

    return False

def futoshiki_solver():
    # Initial grid with 0 representing empty cells
    grid = [
        [0, 3, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 1, 0, 0]
    ]

    # Define the constraints
    constraints = {
        'horizontal': [
            (0, 1, 2, '<'), (0, 2, 3, '>'),
            (1, 0, 1, '>'), (1, 2, 3, '<'),
            (2, 1, 2, '>'), (2, 2, 3, '<')
        ],
        'vertical': [
            (0, 0, 1, '∧'), (1, 0, 1, '∨'), (2, 0, 1, '∨'),
            (2, 1, 2, '∨'), (0, 2, 3, '∧')
        ]
    }

    if solve_futoshiki(grid, 0, 0, constraints):
        return grid
    else:
        return None

solution = futoshiki_solver()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")