def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for j in range(4):
        if grid[row][j] == num:
            return False
    
    # Check column
    for i in range(4):
        if grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:  # Check left constraint
        if h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
            if not (grid[row][col-1] > num):
                return False
        if h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
            if not (grid[row][col-1] < num):
                return False
    
    if col < 3:  # Check right constraint
        if h_constraints[row][col] == '<' and grid[row][col+1] != 0:
            if not (num < grid[row][col+1]):
                return False
        if h_constraints[row][col] == '>' and grid[row][col+1] != 0:
            if not (num > grid[row][col+1]):
                return False
    
    # Check vertical constraints
    if row > 0:  # Check upper constraint
        if v_constraints[row-1][col] == 'v' and grid[row-1][col] != 0:
            if not (grid[row-1][col] > num):
                return False
        if v_constraints[row-1][col] == '^' and grid[row-1][col] != 0:
            if not (grid[row-1][col] < num):
                return False
    
    if row < 3:  # Check lower constraint
        if v_constraints[row][col] == '^' and grid[row+1][col] != 0:
            if not (num < grid[row+1][col]):
                return False
        if v_constraints[row][col] == 'v' and grid[row+1][col] != 0:
            if not (num > grid[row+1][col]):
                return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                for num in range(1, 5):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize the puzzle
grid = [
    [0, 0, 0, 0],
    [0, 3, 0, 4],
    [0, 4, 3, 2],
    [0, 0, 0, 0]
]

# Horizontal constraints (> and <)
h_constraints = [
    ['>', '', '>', ''],
    ['', '', '', ''],
    ['', '', '', ''],
    ['', '', '', '']
]

# Vertical constraints (^ and v)
v_constraints = [
    ['', '', '', ''],
    ['', '', '', ''],
    ['', '', '', '']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    for row in grid:
        print(' '.join(map(str, row)))