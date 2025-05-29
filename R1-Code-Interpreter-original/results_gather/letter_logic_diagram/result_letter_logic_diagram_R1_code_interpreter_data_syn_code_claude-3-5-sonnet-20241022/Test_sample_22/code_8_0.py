def is_valid(grid, row, col, letter, initial_grid):
    # Must match initial grid if position was pre-filled
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

    # Check minor diagonal - must be 'c'
    if row + col == 6 and letter != 'c':
        return False

    return True

def find_next_empty(grid, initial_grid):
    # First try to fill pre-filled positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and initial_grid[i][j] != '':
                return (i, j)
    
    # Then try diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Finally try remaining positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    
    return None

def solve(grid, initial_grid):
    pos = find_next_empty(grid, initial_grid)
    if not pos:
        return True

    row, col = pos

    # If this is a pre-filled position, only try that value
    if initial_grid[row][col] != '':
        letter = initial_grid[row][col]
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
        return False

    # If on minor diagonal, only try 'c'
    if row + col == 6:
        if is_valid(grid, row, col, 'c', initial_grid):
            grid[row][col] = 'c'
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
    ['c','','','','','',''],
    ['','','e','','','',''],
    ['d','e','g','','f','c','b'],
    ['e','g','','f','c','','d'],
    ['','a','f','','','',''],
    ['','f','c','b','','e','g'],
    ['f','c','','','e','g','']
]

# Verify initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            if i + j == 6 and initial_grid[i][j] != 'c':
                valid = False
                break

# Create working copy
grid = [row[:] for row in initial_grid]

if valid and solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")