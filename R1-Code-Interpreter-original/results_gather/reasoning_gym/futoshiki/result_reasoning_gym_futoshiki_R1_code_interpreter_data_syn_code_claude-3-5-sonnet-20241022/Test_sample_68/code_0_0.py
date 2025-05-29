def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(5) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    for (r, c1, c2, greater) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if greater and not (num > grid[r][c2]):
                    return False
                if not greater and not (num < grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if greater and not (grid[r][c1] > num):
                    return False
                if not greater and not (grid[r][c1] < num):
                    return False

    # Check vertical constraints
    for (r1, r2, c, greater) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if greater and not (num > grid[r2][c]):
                    return False
                if not greater and not (num < grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if greater and not (grid[r1][c] > num):
                    return False
                if not greater and not (grid[r1][c] < num):
                    return False
                    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                for num in range(1, 6):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 3, 2, 0],
    [0, 0, 0, 4, 3],
    [0, 0, 0, 1, 0],
    [0, 3, 0, 0, 0],
    [0, 0, 1, 0, 0]
]

# Horizontal constraints: (row, col1, col2, is_greater)
# True means col1 > col2, False means col1 < col2
h_constraints = [
    (1, 3, 4, True),
    (2, 0, 1, True),
    (3, 0, 1, True),
    (4, 3, 4, True)
]

# Vertical constraints: (row1, row2, col, is_greater)
# True means row1 > row2, False means row1 < row2
v_constraints = [
    (0, 1, 0, True),
    (2, 3, 0, True),
    (1, 2, 4, True)
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    for row in grid:
        print(' '.join(map(str, row)))