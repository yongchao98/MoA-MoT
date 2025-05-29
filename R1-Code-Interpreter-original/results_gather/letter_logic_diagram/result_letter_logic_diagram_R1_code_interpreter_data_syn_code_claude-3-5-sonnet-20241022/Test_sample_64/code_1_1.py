def is_valid(grid, row, col, letter, initial_grid):
    # Check initial constraints first
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal - must be 'a'
    if row + col == 6 and letter != 'a':
        return False
    
    return True

def solve(grid, initial_grid):
    # Find first empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, must be 'a'
    if row + col == 6:
        if is_valid(grid, row, col, 'a', initial_grid):
            grid[row][col] = 'a'
            if solve(grid, initial_grid):
                return True
        grid[row][col] = ''
        return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['g', '', '', 'd', '', '', ''],
    ['', 'b', 'd', '', '', '', ''],
    ['', '', '', '', 'e', 'g', 'a'],
    ['', '', 'f', '', '', 'a', ''],
    ['', 'f', 'e', '', '', 'b', ''],
    ['f', '', '', '', 'b', 'd', 'c'],
    ['', 'g', 'a', '', '', '', 'f']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

# Fill all minor diagonal positions with 'a'
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = 'a'

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")