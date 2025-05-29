def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(8) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    for (r, c1, c2, op) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if op == '>' and not (num > grid[r][c2]):
                    return False
                if op == '<' and not (num < grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if op == '>' and not (grid[r][c1] > num):
                    return False
                if op == '<' and not (grid[r][c1] < num):
                    return False
    
    # Check vertical constraints
    for (r1, r2, c, op) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if op == '>' and not (num > grid[r2][c]):
                    return False
                if op == '<' and not (num < grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if op == '>' and not (grid[r1][c] > num):
                    return False
                if op == '<' and not (grid[r1][c] < num):
                    return False
    
    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 9):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty)
grid = [
    [0, 7, 0, 2, 0, 0, 1, 5],
    [7, 0, 0, 0, 0, 4, 0, 2],
    [0, 5, 7, 0, 1, 0, 4, 6],
    [6, 0, 0, 3, 4, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 0],
    [3, 2, 6, 0, 7, 0, 5, 0],
    [0, 4, 2, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 5, 0, 8, 1]
]

# Horizontal constraints: (row, col1, col2, operator)
h_constraints = [
    (0, 2, 3, '>'),
    (3, 4, 5, '<'),
    (4, 0, 1, '<'),
    (5, 2, 3, '>'),
    (5, 4, 5, '>'),
    (6, 0, 1, '>'),
    (7, 4, 5, '>')
]

# Vertical constraints: (row1, row2, col, operator)
v_constraints = [
    (1, 2, 7, '^'),
    (3, 4, 1, 'v'),
    (5, 6, 4, 'v')
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    for row in grid:
        print('   '.join(str(x) for x in row))
else:
    print("No solution exists")