def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col < 5 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    if col < 5 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row < 5 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    if row < 5 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                for num in range(1, 7):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 0, 0, 3, 1],
    [0, 0, 5, 0, 0, 0],
    [0, 0, 2, 0, 0, 0],
    [5, 0, 3, 2, 6, 0],
    [0, 0, 0, 5, 1, 0],
    [4, 0, 0, 0, 0, 0]
]

# Horizontal constraints (> and <)
h_constraints = [
    ['>','','','',''],
    ['','','','',''],
    ['','','','',''],
    ['','','','',''],
    ['','','','',''],
    ['','','','','']
]

# Vertical constraints (∧ and ∨)
v_constraints = [
    ['','','','','',''],
    ['','','','','',''],
    ['','','∧','','',''],
    ['','','','','',''],
    ['','','','','','']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    for row in grid:
        print(' '.join(map(str, row)))