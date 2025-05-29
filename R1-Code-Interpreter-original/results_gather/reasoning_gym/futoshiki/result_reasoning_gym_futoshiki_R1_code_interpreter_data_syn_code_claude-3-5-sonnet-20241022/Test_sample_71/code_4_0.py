def print_grid(grid):
    for row in grid:
        print(' '.join(str(x) if x != 0 else '_' for x in row))

def get_domain(grid, row, col, h_constraints, v_constraints):
    if grid[row][col] != 0:
        return set()
    
    domain = set(range(1, 10))
    
    # Remove numbers used in the same row
    for j in range(9):
        if grid[row][j] != 0:
            domain.discard(grid[row][j])
    
    # Remove numbers used in the same column
    for i in range(9):
        if grid[i][col] != 0:
            domain.discard(grid[i][col])
    
    # Apply horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0:
            domain = {x for x in domain if x > grid[row][col-1]}
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0:
            domain = {x for x in domain if x < grid[row][col-1]}
    if col < 8 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0:
            domain = {x for x in domain if x < grid[row][col+1]}
    if col < 8 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0:
            domain = {x for x in domain if x > grid[row][col+1]}
    
    # Apply vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0:
            domain = {x for x in domain if x > grid[row-1][col]}
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0:
            domain = {x for x in domain if x < grid[row-1][col]}
    if row < 8 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0:
            domain = {x for x in domain if x < grid[row+1][col]}
    if row < 8 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0:
            domain = {x for x in domain if x > grid[row+1][col]}
    
    return domain

def select_unassigned_variable(grid, h_constraints, v_constraints):
    min_domain_size = 10
    best_cell = None
    best_domain = None
    
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                domain = get_domain(grid, i, j, h_constraints, v_constraints)
                if len(domain) < min_domain_size:
                    min_domain_size = len(domain)
                    best_cell = (i, j)
                    best_domain = domain
                    if min_domain_size == 1:
                        return best_cell, best_domain
    
    return best_cell, best_domain

def solve(grid, h_constraints, v_constraints):
    # Find cell with smallest domain
    cell_and_domain = select_unassigned_variable(grid, h_constraints, v_constraints)
    if not cell_and_domain[0]:  # No empty cells left
        return True
    
    row, col = cell_and_domain[0]
    domain = cell_and_domain[1]
    
    if not domain:  # No valid values available
        return False
    
    # Try each value in the domain
    for value in sorted(domain):
        grid[row][col] = value
        if solve(grid, h_constraints, v_constraints):
            return True
        grid[row][col] = 0
    
    return False

# Initial grid with given values
grid = [
    [7,0,8,9,0,6,5,0,0],
    [5,0,0,0,0,0,0,0,0],
    [0,0,0,0,9,7,2,0,3],
    [3,0,1,0,0,4,0,0,0],
    [2,1,0,3,0,0,9,0,6],
    [9,7,6,0,0,0,0,3,0],
    [6,0,5,4,7,3,0,1,0],
    [0,5,0,0,4,0,0,7,0],
    [0,0,0,5,0,0,0,0,8]
]

# Horizontal constraints
h_constraints = [
    [' ','<',' ','>',' ','<',' ','>',' '],
    ['>',' ','<',' ','<',' ',' ',' ',' '],
    ['<',' ',' ',' ',' ',' ',' ',' ',' '],
    [' ',' ','<',' ',' ',' ','<',' ',' '],
    [' ','<',' ',' ',' ','<',' ',' ',' '],
    [' ','>',' ',' ',' ',' ','>',' ',' '],
    [' ',' ',' ',' ','>',' ',' ',' ',' '],
    [' ',' ','>',' ',' ',' ',' ',' ',' '],
    [' ',' ',' ',' ',' ',' ','>',' ',' ']
]

# Vertical constraints
v_constraints = [
    [' ','∨',' ','∧','∧',' ',' ','∨',' '],
    ['∨',' ','∧',' ',' ',' ',' ',' ',' '],
    [' ',' ',' ',' ',' ',' ',' ','∧',' '],
    [' ',' ',' ',' ','∨',' ','∨',' ',' '],
    [' ',' ',' ',' ',' ',' ',' ',' ',' '],
    ['∧',' ',' ','∧',' ','∨',' ','∧','∨'],
    [' ',' ',' ',' ',' ',' ','∧',' ',' '],
    [' ',' ',' ',' ',' ',' ',' ',' ',' '],
    [' ',' ',' ',' ',' ',' ',' ',' ',' ']
]

if solve(grid, h_constraints, v_constraints):
    print_grid(grid)
else:
    print("No solution exists")