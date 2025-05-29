def print_solution(grid):
    for row in grid:
        print(' '.join(map(str, row)))

def valid_move(grid, row, col, num, h_constraints, v_constraints):
    # Check row and column
    for i in range(6):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '<':
            if grid[row][col-1] != 0 and grid[row][col-1] >= num:
                return False
        elif h_constraints[row][col-1] == '>':
            if grid[row][col-1] != 0 and grid[row][col-1] <= num:
                return False
    
    if col < 5:
        if h_constraints[row][col] == '<':
            if grid[row][col+1] != 0 and num >= grid[row][col+1]:
                return False
        elif h_constraints[row][col] == '>':
            if grid[row][col+1] != 0 and num <= grid[row][col+1]:
                return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨':
            if grid[row-1][col] != 0 and grid[row-1][col] <= num:
                return False
        elif v_constraints[row-1][col] == '∧':
            if grid[row-1][col] != 0 and grid[row-1][col] >= num:
                return False
    
    if row < 5:
        if v_constraints[row][col] == '∨':
            if grid[row+1][col] != 0 and num <= grid[row+1][col]:
                return False
        elif v_constraints[row][col] == '∧':
            if grid[row+1][col] != 0 and num >= grid[row+1][col]:
                return False
    
    return True

def find_empty(grid):
    # Find cells with minimum possible values first
    min_possibilities = 10
    best_pos = None
    
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                return (i, j)  # Return first empty cell
    return None

def solve_puzzle(grid, h_constraints, v_constraints):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    for num in range(1, 7):
        if valid_move(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_puzzle(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid with given numbers
initial_grid = [
    [0, 0, 0, 2, 4, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 4, 0, 5, 0, 0],
    [0, 6, 3, 0, 0, 0],
    [0, 0, 6, 0, 3, 5],
    [0, 0, 5, 0, 0, 4]
]

# Horizontal constraints
h_constraints = [
    ['<', ' ', ' ', '<', ' '],
    [' ', ' ', ' ', '>', ' '],
    [' ', ' ', ' ', ' ', ' '],
    [' ', ' ', '>', '<', ' '],
    [' ', ' ', ' ', ' ', ' '],
    ['>', ' ', ' ', ' ', ' ']
]

# Vertical constraints
v_constraints = [
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', '∨', ' ', '∨', ' ', ' ']
]

# Pre-validate initial numbers
valid = True
for i in range(6):
    for j in range(6):
        if initial_grid[i][j] != 0:
            num = initial_grid[i][j]
            initial_grid[i][j] = 0
            if not valid_move(initial_grid, i, j, num, h_constraints, v_constraints):
                valid = False
            initial_grid[i][j] = num

if valid and solve_puzzle(initial_grid, h_constraints, v_constraints):
    print_solution(initial_grid)