def print_grid(grid):
    for row in grid:
        print(" ".join(map(str, row)))

def is_valid_number(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    for i in range(8):
        if grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
            if not (grid[row][col-1] > num):
                return False
        if h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
            if not (grid[row][col-1] < num):
                return False
    
    if col < 7:
        if h_constraints[row][col] == '>' and grid[row][col+1] != 0:
            if not (num > grid[row][col+1]):
                return False
        if h_constraints[row][col] == '<' and grid[row][col+1] != 0:
            if not (num < grid[row][col+1]):
                return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
            if not (grid[row-1][col] > num):
                return False
        if v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
            if not (grid[row-1][col] < num):
                return False
    
    if row < 7:
        if v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
            if not (num > grid[row+1][col]):
                return False
        if v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
            if not (num < grid[row+1][col]):
                return False
    
    return True

def find_empty_location(grid):
    min_possibilities = 9
    best_pos = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)  # Return first empty cell
    return None

def solve(grid, h_constraints, v_constraints):
    empty = find_empty_location(grid)
    
    if not empty:
        return True
    
    row, col = empty
    
    # Special handling for the cell above the '∧' constraint
    if row == 2 and col == 6:  # Position above the ∧ constraint
        # Try smaller numbers first
        numbers = range(1, 5)  # Restrict to smaller numbers
    else:
        numbers = range(1, 9)
    
    for num in numbers:
        if is_valid_number(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            
            if solve(grid, h_constraints, v_constraints):
                return True
            
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 0, 4, 0, 5, 6, 8],
    [6, 0, 0, 0, 0, 3, 5, 2],
    [8, 0, 3, 0, 7, 1, 4, 0],
    [0, 0, 8, 1, 4, 0, 0, 0],
    [0, 0, 1, 0, 5, 0, 0, 7],
    [4, 0, 0, 0, 0, 0, 8, 0],
    [0, 7, 2, 0, 0, 0, 0, 0],
    [0, 2, 0, 8, 0, 0, 0, 0]
]

# Horizontal constraints
h_constraints = [
    ['', '', '', '', '', '', ''],
    ['', '>', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '']
]

# Vertical constraints
v_constraints = [
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '∧', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    print_grid(grid)
else:
    print("No solution exists")