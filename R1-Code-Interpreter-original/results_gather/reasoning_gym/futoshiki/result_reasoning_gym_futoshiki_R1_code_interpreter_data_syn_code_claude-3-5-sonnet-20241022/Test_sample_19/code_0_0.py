def is_valid(grid, row, col, num, constraints_h, constraints_v):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if constraints_h[row][col-1] == '>' and num >= grid[row][col-1]:
            return False
        if constraints_h[row][col-1] == '<' and num <= grid[row][col-1]:
            return False
    if col < 4:
        if constraints_h[row][col] == '<' and num >= grid[row][col+1]:
            return False
        if constraints_h[row][col] == '>' and num <= grid[row][col+1]:
            return False
    
    # Check vertical constraints
    if row > 0:
        if constraints_v[row-1][col] == '∧' and num >= grid[row-1][col]:
            return False
        if constraints_v[row-1][col] == '∨' and num <= grid[row-1][col]:
            return False
    if row < 4:
        if constraints_v[row][col] == '∨' and num >= grid[row+1][col]:
            return False
        if constraints_v[row][col] == '∧' and num <= grid[row+1][col]:
            return False
    
    return True

def find_empty(grid):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid, constraints_h, constraints_v):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 6):
        if is_valid(grid, row, col, num, constraints_h, constraints_v):
            grid[row][col] = num
            if solve_futoshiki(grid, constraints_h, constraints_v):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the grid (0 represents empty cells)
grid = [
    [0, 0, 0, 5, 0],
    [0, 0, 0, 0, 4],
    [4, 0, 0, 0, 0],
    [0, 0, 2, 0, 0],
    [5, 0, 4, 0, 0]
]

# Initialize horizontal constraints ('' means no constraint)
constraints_h = [
    ['', '', '', '', ''],
    ['>', '', '>', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '>', '', '']
]

# Initialize vertical constraints
constraints_v = [
    ['', '', '', '', ''],
    ['', '', '∧', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', '∧'],
    ['', '', '', '', '']
]

if solve_futoshiki(grid, constraints_h, constraints_v):
    result = ""
    for i in range(5):
        row = ""
        for j in range(5):
            row += str(grid[i][j])
            if j < 4:
                row += " " + constraints_h[i][j] + " " if constraints_h[i][j] else "   "
        result += row + "\n"
        if i < 4:
            for j in range(5):
                result += constraints_v[i][j] + "   " if constraints_v[i][j] else "    "
            result += "\n"
    print(result)
else:
    print("No solution exists")