def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(9) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    for (r, c1, c2, type) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if type == '<' and not (num < grid[r][c2]):
                    return False
                if type == '>' and not (num > grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if type == '<' and not (grid[r][c1] < num):
                    return False
                if type == '>' and not (grid[r][c1] > num):
                    return False

    # Check vertical constraints
    for (r1, r2, c, type) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if type == '∨' and not (num > grid[r2][c]):
                    return False
                if type == '∧' and not (num < grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if type == '∨' and not (grid[r1][c] > num):
                    return False
                if type == '∧' and not (grid[r1][c] < num):
                    return False
    
    return True

def find_empty(grid):
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 10):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [8, 5, 9, 0, 4, 0, 1, 0, 0],
    [0, 4, 2, 5, 0, 0, 3, 0, 0],
    [7, 1, 0, 0, 0, 9, 0, 4, 0],
    [0, 0, 0, 0, 0, 6, 7, 3, 0],
    [3, 7, 0, 0, 9, 1, 0, 0, 0],
    [0, 0, 0, 2, 8, 0, 0, 0, 0],
    [6, 0, 7, 0, 0, 0, 0, 5, 0],
    [0, 9, 6, 0, 2, 0, 0, 0, 8],
    [0, 0, 5, 7, 0, 2, 0, 0, 0]
]

# Horizontal constraints (row, col1, col2, type)
h_constraints = [
    (0, 2, 3, '<'), (0, 4, 5, '>'), (0, 5, 6, '>'), (0, 6, 7, '>'),
    (1, 0, 1, '<'), (1, 6, 7, '>'), (1, 7, 8, '>'),
    (2, 3, 4, '<'), (2, 5, 6, '>'), (2, 6, 7, '>'),
    (3, 1, 2, '>'), (3, 3, 4, '>'),
    (4, 2, 3, '>'), (4, 3, 4, '<'),
    (6, 4, 5, '<'), (6, 6, 7, '>'),
    (7, 2, 3, '>'), (7, 3, 4, '>'), (7, 6, 7, '>'),
]

# Vertical constraints (row1, row2, col, type)
v_constraints = [
    (1, 2, 4, '∨'), (3, 4, 3, '∨'), (3, 4, 6, '∨'),
    (4, 5, 6, '∧'), (6, 7, 0, '∨'), (7, 8, 1, '∨'),
    (3, 4, 8, '∧')
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    for row in grid:
        print(row)
else:
    print("No solution exists")