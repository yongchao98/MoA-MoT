def print_grid(grid):
    for row in grid:
        print("   ".join(map(str, row)))

def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for j in range(8):
        if grid[row][j] == num:
            return False
            
    # Check column
    for i in range(8):
        if grid[i][col] == num:
            return False
            
    # Check horizontal constraints
    # Left constraint
    if col > 0:
        if h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
            if not (grid[row][col-1] > num):
                return False
        elif h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
            if not (grid[row][col-1] < num):
                return False
                
    # Right constraint
    if col < 7:
        if h_constraints[row][col] == '>' and grid[row][col+1] != 0:
            if not (num > grid[row][col+1]):
                return False
        elif h_constraints[row][col] == '<' and grid[row][col+1] != 0:
            if not (num < grid[row][col+1]):
                return False
                
    # Check vertical constraints
    # Upper constraint
    if row > 0:
        if v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
            if not (grid[row-1][col] > num):
                return False
        elif v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
            if not (grid[row-1][col] < num):
                return False
                
    # Lower constraint
    if row < 7:
        if v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
            if not (num > grid[row+1][col]):
                return False
        elif v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
            if not (num < grid[row+1][col]):
                return False
                
    return True

def find_empty(grid):
    min_options = 9
    best_pos = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    
    if not empty:
        return True
        
    row, col = empty
    
    # Try values in different orders based on position
    values = list(range(1, 9))
    if row < 4:
        values.reverse()
    
    for num in values:
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            
            if solve(grid, h_constraints, v_constraints):
                return True
                
            grid[row][col] = 0
            
    return False

# Initialize puzzle
grid = [
    [2,0,0,4,0,5,6,0],
    [0,8,2,0,0,0,0,0],
    [8,0,4,7,0,3,6,2],
    [0,1,0,0,0,0,0,7],
    [0,7,5,0,0,0,0,0],
    [1,0,6,0,0,0,7,0],
    [0,4,7,3,0,5,0,0],
    [0,0,0,0,0,0,1,5]
]

h_constraints = [
    ['','>','<','','','>','',''],
    ['','','>','','','','',''],
    ['','','','','','','',''],
    ['>','','','<','>','','',''],
    ['','','<','','','<','',''],
    ['<','','','','','','<',''],
    ['','','','','','','',''],
    ['','>','','','','','','']
]

v_constraints = [
    ['','∨','','','∨','','∨',''],
    ['','','','','','','',''],
    ['∨','','','','','','',''],
    ['∨','','','∨','','','∨',''],
    ['∨','','','∧','','','',''],
    ['','∧','∧','','∧','∧','',''],
    ['∨','∧','','','','','∨',''],
    ['','','','','','','','']
]

# Verify initial state
valid_initial = True
for i in range(8):
    for j in range(8):
        if grid[i][j] != 0:
            val = grid[i][j]
            grid[i][j] = 0
            if not is_valid(grid, i, j, val, h_constraints, v_constraints):
                valid_initial = False
            grid[i][j] = val

if not valid_initial:
    print("Invalid initial state")
else:
    if solve(grid, h_constraints, v_constraints):
        for i in range(8):
            row = []
            for j in range(8):
                row.append(str(grid[i][j]))
            print("   ".join(row))
    else:
        print("No solution exists")