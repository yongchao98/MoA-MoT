def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(8) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
        if not (grid[row][col-1] < num):
            return False
    if col > 0 and h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
        if not (grid[row][col-1] > num):
            return False
    if col < 7 and h_constraints[row][col] == '<' and grid[row][col+1] != 0:
        if not (num < grid[row][col+1]):
            return False
    if col < 7 and h_constraints[row][col] == '>' and grid[row][col+1] != 0:
        if not (num > grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
        if not (grid[row-1][col] < num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
        if not (grid[row-1][col] > num):
            return False
    if row < 7 and v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
        if not (num < grid[row+1][col]):
            return False
    if row < 7 and v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
        if not (num > grid[row+1][col]):
            return False
    
    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 9):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the puzzle
grid = [
    [3, 5, 0, 0, 0, 7, 0, 0],
    [0, 2, 5, 0, 6, 0, 0, 7],
    [0, 0, 0, 6, 1, 0, 7, 0],
    [0, 0, 0, 2, 0, 0, 8, 0],
    [5, 7, 0, 0, 0, 2, 0, 0],
    [0, 0, 6, 0, 4, 0, 0, 2],
    [8, 0, 0, 0, 0, 0, 2, 3],
    [0, 0, 0, 0, 3, 6, 4, 0]
]

# Horizontal constraints ('' means no constraint)
h_constraints = [
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '<', '>', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '<', '', '', '', ''],
    ['<', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '<', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '>', '']
]

# Vertical constraints
v_constraints = [
    ['', '', '', '', '', '', '', ''],
    ['', '∧', '', '', '', '', '', ''],
    ['', '', '∧', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '∧', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j])
            if j < 7:
                row += "   " if h_constraints[i][j] == '' else f" {h_constraints[i][j]} "
        result += row + "\n"
        if i < 7:
            v_row = ""
            for j in range(8):
                v_row += v_constraints[i][j] if v_constraints[i][j] != '' else " "
                v_row += "   "
            if v_row.strip():
                result += v_row + "\n"
    print(result)