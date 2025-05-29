def print_grid(grid):
    for i in range(9):
        row = ""
        for j in range(9):
            row += str(grid[i][j]) if grid[i][j] != 0 else "_"
            row += "   "
        print(row.rstrip())

def valid_in_row_col(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    # Check column
    if num in [grid[i][col] for i in range(9)]:
        return False
    return True

def check_horizontal_constraint(grid, row, col1, col2, type, num=None):
    val1 = num if num is not None else grid[row][col1]
    val2 = grid[row][col2]
    
    if val1 == 0 or val2 == 0:
        return True
        
    if type == '<':
        return val1 < val2
    return val1 > val2

def check_vertical_constraint(grid, row1, row2, col, type, num=None):
    val1 = num if num is not None else grid[row1][col]
    val2 = grid[row2][col]
    
    if val1 == 0 or val2 == 0:
        return True
        
    if type == '∨':
        return val1 > val2
    return val1 < val2

def is_valid(grid, row, col, num, h_constraints, v_constraints):
    if not valid_in_row_col(grid, row, col, num):
        return False

    # Check horizontal constraints
    for r, c1, c2, type in h_constraints:
        if r == row:
            if c1 == col:
                if not check_horizontal_constraint(grid, r, col, c2, type, num):
                    return False
            elif c2 == col:
                if not check_horizontal_constraint(grid, r, c1, col, type, grid[r][c1]):
                    return False

    # Check vertical constraints
    for r1, r2, c, type in v_constraints:
        if c == col:
            if r1 == row:
                if not check_vertical_constraint(grid, row, r2, c, type, num):
                    return False
            elif r2 == row:
                if not check_vertical_constraint(grid, r1, row, c, type, grid[r1][c]):
                    return False

    return True

def find_empty(grid):
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 10):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            
            if solve(grid, h_constraints, v_constraints):
                return True
                
            grid[row][col] = 0
    
    return False

# Initial grid
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

# Horizontal constraints
h_constraints = [
    (0, 2, 3, '<'), (0, 4, 5, '>'), (0, 5, 6, '>'), (0, 6, 7, '>'),
    (1, 0, 1, '<'), (1, 6, 7, '>'), (1, 7, 8, '>'),
    (2, 3, 4, '<'), (2, 5, 6, '>'), (2, 6, 7, '>'),
    (3, 1, 2, '>'), (3, 3, 4, '>'),
    (4, 2, 3, '>'), (4, 3, 4, '<'),
    (6, 4, 5, '<'), (6, 6, 7, '>'),
    (7, 2, 3, '>'), (7, 3, 4, '>'), (7, 6, 7, '>')
]

# Vertical constraints
v_constraints = [
    (1, 2, 4, '∨'),
    (3, 4, 3, '∨'),
    (3, 4, 6, '∨'),
    (4, 5, 6, '∧'),
    (6, 7, 0, '∨'),
    (7, 8, 1, '∨'),
    (3, 4, 8, '∧')
]

print("Initial grid:")
print_grid(grid)
print("\nSolving...")
if solve(grid, h_constraints, v_constraints):
    print("\nSolution found:")
    print_grid(grid)
else:
    print("No solution exists")