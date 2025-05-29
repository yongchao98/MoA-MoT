def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row uniqueness
    for j in range(5):
        if grid[row][j] == num:
            return False
    
    # Check column uniqueness
    for i in range(5):
        if grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    # Left constraint
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
            
    # Right constraint
    if col < 4 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
    if col < 4 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    
    # Check vertical constraints
    # Upper constraint
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
            
    # Lower constraint
    if row < 4 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    if row < 4 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    
    return True

def find_empty(grid):
    # Find cells that need to be filled, prioritizing cells with constraints
    min_options = 6
    best_cell = None
    
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                valid_count = 0
                for num in range(1, 6):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        valid_count += 1
                if valid_count < min_options:
                    min_options = valid_count
                    best_cell = (i, j)
    
    return best_cell

def solve(grid, h_constraints, v_constraints):
    cell = find_empty(grid)
    if not cell:
        return True
    
    row, col = cell
    # Try numbers in an optimized order based on constraints
    numbers = list(range(1, 6))
    if col < 4 and h_constraints[row][col] == '<':
        numbers.sort()  # Try smaller numbers first
    elif col > 0 and h_constraints[row][col-1] == '<':
        numbers.sort(reverse=True)  # Try larger numbers first
    elif row < 4 and v_constraints[row][col] == '∨':
        numbers.sort(reverse=True)  # Try larger numbers first
    
    for num in numbers:
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 4, 5, 3],
    [0, 4, 0, 0, 0],
    [0, 0, 5, 4, 2],
    [0, 0, 2, 0, 0],
    [0, 0, 0, 0, 0]
]

# Horizontal constraints
h_constraints = [
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '<', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '<', '']
]

# Vertical constraints
v_constraints = [
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['∨', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    result = ""
    for i in range(5):
        result += "  ".join(str(x) for x in grid[i]) + "\n"
    print(result)