def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col < 4 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    if col < 4 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row < 4 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    if row < 4 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    
    return True

def solve(grid, h_constraints, v_constraints):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                for num in range(1, 6):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize the grid
grid = [
    [0, 0, 0, 1, 4],
    [0, 0, 3, 0, 0],
    [0, 0, 1, 0, 3],
    [0, 0, 0, 0, 2],
    [0, 0, 0, 0, 0]
]

# Initialize horizontal constraints (> and <)
h_constraints = [
    ['', '', '', '', ''],
    ['>', '', '', '', ''],
    ['', '>', '', '', ''],
    ['', '', '', '>', ''],
    ['', '>', '', '', '']
]

# Initialize vertical constraints (∧ and ∨)
v_constraints = [
    ['', '', '', '∧', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['∨', '', '', '∨', ''],
    ['', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    result = ""
    for i in range(5):
        row = ""
        for j in range(5):
            row += str(grid[i][j])
            if j < 4:
                row += "   " + (h_constraints[i][j] if h_constraints[i][j] else " ")  + "   "
        result += row + "\n"
        if i < 4:
            for j in range(5):
                result += (v_constraints[i][j] if v_constraints[i][j] else " ") + "    "
            result += "\n"
    print(result)