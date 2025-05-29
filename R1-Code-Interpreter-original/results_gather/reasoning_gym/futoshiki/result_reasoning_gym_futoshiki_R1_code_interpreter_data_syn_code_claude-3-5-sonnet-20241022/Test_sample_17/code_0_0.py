def is_valid(grid, row, col, num, constraints_h, constraints_v):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and constraints_h[row][col-1] == '<' and grid[row][col-1] >= num:
        return False
    if col > 0 and constraints_h[row][col-1] == '>' and grid[row][col-1] <= num:
        return False
    if col < 3 and constraints_h[row][col] == '<' and num >= grid[row][col+1]:
        return False
    if col < 3 and constraints_h[row][col] == '>' and num <= grid[row][col+1]:
        return False
    
    # Check vertical constraints
    if row > 0 and constraints_v[row-1][col] == '∧' and grid[row-1][col] <= num:
        return False
    if row > 0 and constraints_v[row-1][col] == '∨' and grid[row-1][col] >= num:
        return False
    if row < 3 and constraints_v[row][col] == '∧' and num <= grid[row+1][col]:
        return False
    if row < 3 and constraints_v[row][col] == '∨' and num >= grid[row+1][col]:
        return False
    
    return True

def find_empty(grid):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid, constraints_h, constraints_v):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 5):
        if is_valid(grid, row, col, num, constraints_h, constraints_v):
            grid[row][col] = num
            if solve_futoshiki(grid, constraints_h, constraints_v):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the grid
grid = [
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [4, 1, 0, 0],
    [0, 0, 4, 0]
]

# Initialize constraints ('' means no constraint)
constraints_h = [
    ['', '', ''],
    ['', '', ''],
    ['', '', '<'],
    ['', '', '']
]

constraints_v = [
    ['', '∧', '', ''],
    ['', '', '', ''],
    ['', '', '', '']
]

if solve_futoshiki(grid, constraints_h, constraints_v):
    result = ""
    for i in range(4):
        result += "  ".join(str(x) for x in grid[i]) + "\n"
        if i < 3:
            result += "    " + "  ".join(constraints_v[i]) + "\n" if any(constraints_v[i]) else "\n"
    print(result.rstrip())
else:
    print("No solution exists")