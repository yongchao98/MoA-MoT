def find_empty(grid):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                return (i, j)
    return None

def is_valid_move(grid, row, col, num, constraints_h, constraints_v):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    for i in range(5):
        if grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and constraints_h[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col > 0 and constraints_h[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col < 4 and constraints_h[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
    if col < 4 and constraints_h[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
            
    # Check vertical constraints
    if row > 0 and constraints_v[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row > 0 and constraints_v[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row < 4 and constraints_v[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    if row < 4 and constraints_v[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    
    return True

def solve(grid, constraints_h, constraints_v):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    for num in range(1, 6):
        if is_valid_move(grid, row, col, num, constraints_h, constraints_v):
            grid[row][col] = num
            if solve(grid, constraints_h, constraints_v):
                return True
            grid[row][col] = 0
    
    return False

# Initialize puzzle
grid = [
    [0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 2, 4, 0],
    [2, 0, 1, 3, 0],
    [5, 0, 0, 0, 0]
]

constraints_h = [
    ['', '', '', '<', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['>', '', '', '<', '']
]

constraints_v = [
    ['∧', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '∨', '', '', '']
]

if solve(grid, constraints_h, constraints_v):
    result = ""
    for i in range(5):
        row = ""
        for j in range(5):
            row += str(grid[i][j])
            if j < 4:
                row += " " + constraints_h[i][j] + " " if constraints_h[i][j] else "   "
        print(row)
        if i < 4:
            for j in range(5):
                print(constraints_v[i][j] + "   " if constraints_v[i][j] else "    ", end="")
            print()
else:
    print("No solution exists")