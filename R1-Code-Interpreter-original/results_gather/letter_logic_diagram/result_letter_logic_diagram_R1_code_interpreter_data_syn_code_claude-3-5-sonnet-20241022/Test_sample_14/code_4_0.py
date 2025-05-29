def print_grid(grid):
    for row in grid:
        print(','.join(row))

def verify_initial_constraints(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def verify_constraints(grid, row, col, letter, initial_grid):
    # Check if this violates initial constraints
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
            
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid, initial_grid):
    # First try to fill cells that are specified in initial_grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and initial_grid[i][j] != '':
                return (i, j)
    
    # Then fill remaining cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial_grid):
    pos = find_empty(grid, initial_grid)
    if not pos:
        return True
    
    row, col = pos
    
    # If this position is specified in initial_grid, only try that letter
    if initial_grid[row][col] != '':
        letter = initial_grid[row][col]
        if verify_constraints(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
        return False
    
    # If this is on minor diagonal, try to use existing minor diagonal letter
    if row + col == 6:
        minor_diag_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                minor_diag_letter = grid[i][j]
                break
        if minor_diag_letter:
            if verify_constraints(grid, row, col, minor_diag_letter, initial_grid):
                grid[row][col] = minor_diag_letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False

    # Try each possible letter
    for letter in 'abcdefg':
        if verify_constraints(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', 'f', '', 'c', '', '', ''],
    ['', 'd', '', 'e', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['c', 'e', 'g', '', 'b', 'f', ''],
    ['', '', '', '', 'f', 'd', ''],
    ['', 'a', 'b', '', '', 'c', ''],
    ['', 'b', 'f', 'd', '', '', 'g']
]

grid = [row[:] for row in initial_grid]

if solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")