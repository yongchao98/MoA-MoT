def print_solution(grid):
    for row in grid:
        print(' '.join(str(x) for x in row))

def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(6)]:
        return False
    
    # Check horizontal constraints
    if col > 0:  # Check left constraint
        if h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
            if not (num > grid[row][col-1]):
                return False
        elif h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
            if not (num < grid[row][col-1]):
                return False
    if col < 5:  # Check right constraint
        if h_constraints[row][col] == '<' and grid[row][col+1] != 0:
            if not (num < grid[row][col+1]):
                return False
        elif h_constraints[row][col] == '>' and grid[row][col+1] != 0:
            if not (num > grid[row][col+1]):
                return False
    
    # Check vertical constraints
    if row > 0:  # Check upper constraint
        if v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
            if not (num < grid[row-1][col]):
                return False
        elif v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
            if not (num > grid[row-1][col]):
                return False
    if row < 5:  # Check lower constraint
        if v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
            if not (num > grid[row+1][col]):
                return False
        elif v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
            if not (num < grid[row+1][col]):
                return False
    
    return True

def solve(grid, h_constraints, v_constraints):
    for row in range(6):
        for col in range(6):
            if grid[row][col] == 0:
                for num in range(1, 7):
                    if is_valid(grid, row, col, num, h_constraints, v_constraints):
                        grid[row][col] = num
                        if solve(grid, h_constraints, v_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [0, 3, 6, 0, 4, 0],
    [0, 6, 2, 4, 0, 0],
    [6, 2, 0, 0, 1, 4],
    [0, 0, 0, 1, 0, 0],
    [0, 4, 0, 3, 2, 0],
    [0, 0, 0, 0, 0, 0]
]

# Horizontal constraints (empty string means no constraint)
h_constraints = [
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', '']
]

# Vertical constraints (empty string means no constraint)
v_constraints = [
    ['', '', '', '', '', ''],
    ['', '', '', '', '', ''],
    ['', '', '', '', '', ''],
    ['', '', '', '', '', ''],
    ['', '', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    print_solution(grid)
else:
    print("No solution exists")