def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Skip checking against itself
    temp = grid[row][col]
    grid[row][col] = 0
    
    # Check row
    for j in range(8):
        if grid[row][j] == num:
            grid[row][col] = temp
            return False
            
    # Check column
    for i in range(8):
        if grid[i][col] == num:
            grid[row][col] = temp
            return False
            
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
            if not (grid[row][col-1] > num):
                grid[row][col] = temp
                return False
        if h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
            if not (grid[row][col-1] < num):
                grid[row][col] = temp
                return False
    
    if col < 7:
        if h_constraints[row][col] == '>' and grid[row][col+1] != 0:
            if not (num > grid[row][col+1]):
                grid[row][col] = temp
                return False
        if h_constraints[row][col] == '<' and grid[row][col+1] != 0:
            if not (num < grid[row][col+1]):
                grid[row][col] = temp
                return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
            if not (grid[row-1][col] > num):
                grid[row][col] = temp
                return False
        if v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
            if not (grid[row-1][col] < num):
                grid[row][col] = temp
                return False
    
    if row < 7:
        if v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
            if not (num > grid[row+1][col]):
                grid[row][col] = temp
                return False
        if v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
            if not (num < grid[row+1][col]):
                grid[row][col] = temp
                return False
    
    grid[row][col] = temp
    return True

def solve(grid, h_constraints, v_constraints):
    empty = None
    # Find first empty cell
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 9):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the puzzle
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

# First verify all initial numbers are valid
valid = True
for i in range(8):
    for j in range(8):
        if grid[i][j] != 0:
            num = grid[i][j]
            if not is_valid(grid, i, j, num, h_constraints, v_constraints):
                valid = False
                break
    if not valid:
        break

if not valid:
    print("Initial configuration is invalid")
else:
    if solve(grid, h_constraints, v_constraints):
        for i in range(8):
            row = []
            for j in range(8):
                row.append(str(grid[i][j]))
            print("   ".join(row))
    else:
        print("No solution exists")