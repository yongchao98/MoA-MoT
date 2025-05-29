def print_solution(grid):
    for row in grid:
        print(' '.join(map(str, row)))

def valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row and column uniqueness
    for x in range(6):
        if grid[row][x] == num or grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '<':
            if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
                return False
        elif h_constraints[row][col-1] == '>':
            if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
                return False
    if col < 5:
        if h_constraints[row][col] == '<':
            if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
                return False
        elif h_constraints[row][col] == '>':
            if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
                return False

    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row < 5 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False

    return True

def find_next_empty(grid):
    # First, try cells with constraints
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid, h_constraints, v_constraints):
    pos = find_next_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    # Try values in different orders based on position
    values = list(range(1, 7))
    if row == 0 and col == 0:  # First position should try smaller numbers first
        values = [1, 2, 3, 4, 5, 6]
    
    for num in values:
        if valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    return False

# Initial grid with given numbers
grid = [
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

if solve(grid, h_constraints, v_constraints):
    print_solution(grid)