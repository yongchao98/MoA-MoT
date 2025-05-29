def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(9):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(9):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col < 8 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    if col < 8 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
            
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row < 8 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    if row < 8 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
            
    return True

def solve(grid, h_constraints, v_constraints):
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                for num in range(1, 10):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize the puzzle
grid = [
    [0, 9, 0, 0, 8, 0, 0, 3, 0],
    [0, 2, 0, 0, 9, 0, 0, 7, 0],
    [2, 0, 0, 8, 0, 0, 0, 4, 3],
    [4, 1, 0, 0, 0, 8, 0, 0, 0],
    [8, 4, 0, 5, 0, 0, 6, 2, 7],
    [0, 7, 0, 1, 5, 6, 0, 0, 4],
    [0, 0, 0, 0, 1, 7, 9, 8, 0],
    [6, 0, 0, 0, 0, 3, 4, 0, 0],
    [9, 6, 0, 7, 0, 2, 8, 5, 0]
]

# Initialize constraints (empty string means no constraint)
h_constraints = [
    ['', '', '', '', '', '', '', '>', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '>', ''],
    ['', '', '', '', '<', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '<', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', '']
]

v_constraints = [
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '∨', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    # Format the output
    result = ""
    for i in range(9):
        row = ""
        for j in range(9):
            row += str(grid[i][j])
            if j < 8:
                if h_constraints[i][j] != '':
                    row += " " + h_constraints[i][j] + " "
                else:
                    row += "   "
        result += row + "\n"
        if i < 8:
            for j in range(9):
                if v_constraints[i][j] != '':
                    result += v_constraints[i][j] + "   "
                else:
                    result += "    "
            result += "\n"
    print(result)