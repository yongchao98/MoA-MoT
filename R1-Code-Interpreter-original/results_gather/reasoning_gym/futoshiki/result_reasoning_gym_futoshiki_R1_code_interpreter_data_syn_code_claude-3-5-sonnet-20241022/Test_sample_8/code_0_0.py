def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '<' and not (grid[row][col-1] < num):
            return False
        if h_constraints[row][col-1] == '>' and not (grid[row][col-1] > num):
            return False
    if col < 6:
        if h_constraints[row][col] == '<' and not (num < grid[row][col+1] if grid[row][col+1] != 0 else True):
            return False
        if h_constraints[row][col] == '>' and not (num > grid[row][col+1] if grid[row][col+1] != 0 else True):
            return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == 'v' and not (grid[row-1][col] < num):
            return False
        if v_constraints[row-1][col] == '^' and not (grid[row-1][col] > num):
            return False
    if row < 6:
        if v_constraints[row][col] == 'v' and not (num < grid[row+1][col] if grid[row+1][col] != 0 else True):
            return False
        if v_constraints[row][col] == '^' and not (num > grid[row+1][col] if grid[row+1][col] != 0 else True):
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return i, j
    return None

def solve(grid, h_constraints, v_constraints):
    find = find_empty(grid)
    if not find:
        return True
    
    row, col = find
    for num in range(1, 8):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 6, 4, 0, 0, 0, 2],
    [7, 0, 0, 0, 1, 0, 0],
    [0, 2, 0, 0, 4, 6, 0],
    [4, 0, 2, 0, 0, 1, 0],
    [2, 0, 6, 0, 3, 0, 0],
    [0, 1, 3, 7, 0, 0, 0],
    [0, 0, 0, 3, 6, 0, 0]
]

# Horizontal constraints ('<' or '>' or '')
h_constraints = [
    ['', '', '<', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '']
]

# Vertical constraints ('v' for down, '^' for up, or '')
v_constraints = [
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'v', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    result = ""
    for i in range(7):
        if i > 0:
            result += "\n"
        for j in range(7):
            result += str(grid[i][j])
            if j < 6:
                result += "   " + (h_constraints[i][j] if h_constraints[i][j] else " ")
        if i < 6 and any(v_constraints[i]):
            result += "\n"
            for j in range(7):
                if v_constraints[i][j]:
                    result += "    " + v_constraints[i][j]
                else:
                    result += "     "
    print(result)