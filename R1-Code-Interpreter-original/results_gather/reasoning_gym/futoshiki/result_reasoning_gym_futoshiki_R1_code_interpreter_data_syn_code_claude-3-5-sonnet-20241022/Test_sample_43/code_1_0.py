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
        if h_constraints[row][col-1] == '>' and num <= grid[row][col-1]:
            return False
        if h_constraints[row][col-1] == '<' and num >= grid[row][col-1]:
            return False
    if col < 6:
        if h_constraints[row][col] == '>' and num <= grid[row][col+1]:
            return False
        if h_constraints[row][col] == '<' and num >= grid[row][col+1]:
            return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨' and num <= grid[row-1][col]:
            return False
        if v_constraints[row-1][col] == '∧' and num >= grid[row-1][col]:
            return False
    if row < 6:
        if v_constraints[row][col] == '∨' and num >= grid[row+1][col]:
            return False
        if v_constraints[row][col] == '∧' and num <= grid[row+1][col]:
            return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                for num in range(1, 8):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize the grid (0 represents empty cells)
grid = [
    [0, 0, 3, 0, 5, 0, 0],
    [0, 5, 0, 3, 0, 0, 1],
    [0, 0, 0, 2, 0, 4, 0],
    [2, 0, 6, 0, 1, 7, 0],
    [5, 0, 1, 6, 0, 0, 7],
    [0, 0, 7, 0, 0, 0, 2],
    [0, 1, 0, 0, 0, 2, 0]
]

# Horizontal constraints ('' means no constraint)
h_constraints = [
    ['>', '', '', '', '', ''],
    ['<', '', '', '', '', ''],
    ['', '', '', '', '', ''],
    ['', '', '>', '>', '', ''],
    ['', '>', '', '', '', '<'],
    ['<', '', '', '', '', ''],
    ['', '', '', '', '', '']
]

# Vertical constraints
v_constraints = [
    ['', '', '', '', '', '', ''],
    ['', '∨', '', '', '', '', ''],
    ['', '', '', '', '', '∧', ''],
    ['', '', '∧', '', '∧', '', '∧'],
    ['∧', '', '', '', '∧', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    for row in grid:
        print(row)