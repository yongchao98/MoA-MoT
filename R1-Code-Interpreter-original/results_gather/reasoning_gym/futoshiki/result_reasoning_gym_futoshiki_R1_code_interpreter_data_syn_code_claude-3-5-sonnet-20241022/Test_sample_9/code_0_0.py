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
    if col > 0 and h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
        if not (grid[row][col-1] < num):
            return False
    if col > 0 and h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
        if not (grid[row][col-1] > num):
            return False
    if col < 6 and h_constraints[row][col] == '<' and grid[row][col+1] != 0:
        if not (num < grid[row][col+1]):
            return False
    if col < 6 and h_constraints[row][col] == '>' and grid[row][col+1] != 0:
        if not (num > grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
        if not (grid[row-1][col] < num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
        if not (grid[row-1][col] > num):
            return False
    if row < 6 and v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
        if not (num < grid[row+1][col]):
            return False
    if row < 6 and v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
        if not (num > grid[row+1][col]):
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 8):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 5, 0, 4, 0, 0],
    [0, 0, 0, 0, 0, 6, 4],
    [7, 0, 0, 0, 0, 4, 0],
    [0, 2, 0, 0, 7, 1, 0],
    [0, 0, 0, 0, 3, 0, 0],
    [0, 1, 3, 2, 0, 0, 0],
    [6, 4, 0, 0, 1, 0, 0]
]

# Horizontal constraints ('' means no constraint)
h_constraints = [
    ['', '', '', '<', '', ''],
    ['<', '', '', '<', '', ''],
    ['', '', '', '', '<', ''],
    ['', '', '', '', '', ''],
    ['', '', '>', '', '', ''],
    ['', '', '', '', '', ''],
    ['', '', '<', '>', '', '']
]

# Vertical constraints
v_constraints = [
    ['∧', '', '∧', '∨', '∧', ''],
    ['∧', '', '', '∧', '', ''],
    ['', '', '∧', '', '', '∨'],
    ['', '', '', '', '', ''],
    ['', '', '', '', '', ''],
    ['∧', '', '', '', '', '']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j])
            if j < 6:
                row += " " + (h_constraints[i][j] if h_constraints[i][j] else " ") + " "
        result += row + "\n"
        if i < 6:
            for j in range(7):
                result += (v_constraints[i][j] if v_constraints[i][j] else " ") + "   "
            result += "\n"
    print(result)