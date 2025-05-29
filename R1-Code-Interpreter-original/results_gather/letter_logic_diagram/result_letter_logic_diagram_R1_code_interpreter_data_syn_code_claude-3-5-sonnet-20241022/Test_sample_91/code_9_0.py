def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_minor_diagonal_options(grid):
    # Get possible letters for minor diagonal based on initial constraints
    options = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            return {grid[i][6-i]}
    return options

def get_available_letters(grid, row, col, used_in_diag):
    used = set()
    # Check row
    for j in range(7):
        if grid[row][j] != '':
            used.add(grid[row][j])
    
    # Check column
    for i in range(7):
        if grid[i][col] != '':
            used.add(grid[i][col])
    
    # If on diagonal, must use diagonal letter
    if row + col == 6:
        return {used_in_diag} - used
    
    # Can't use diagonal letter in other positions
    available = set('abcdefg') - used - {used_in_diag}
    return available

def is_valid_partial(grid):
    # Check rows for duplicates
    for row in grid:
        seen = set()
        for cell in row:
            if cell != '' and cell in seen:
                return False
            seen.add(cell)
    
    # Check columns for duplicates
    for j in range(7):
        seen = set()
        for i in range(7):
            if grid[i][j] != '' and grid[i][j] in seen:
                return False
            seen.add(grid[i][j])
    
    # Check minor diagonal
    diag = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diag.add(grid[i][6-i])
    if len(diag) > 1:
        return False
    
    # Check initial constraints
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    return True

def solve(grid, pos=0, diag_letter=None):
    if pos == 49:  # 7x7 = 49 cells
        return True
    
    row = pos // 7
    col = pos % 7
    
    # Skip pre-filled cells
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, pos + 1, diag_letter)
    
    # If first unfilled diagonal cell, try each possible diagonal letter
    if row + col == 6 and diag_letter is None:
        for letter in get_minor_diagonal_options(grid):
            grid[row][col] = letter
            if is_valid_partial(grid) and solve(grid, pos + 1, letter):
                return True
            grid[row][col] = ''
        return False
    
    # For other cells
    available = get_available_letters(grid, row, col, diag_letter)
    for letter in available:
        grid[row][col] = letter
        if is_valid_partial(grid) and solve(grid, pos + 1, diag_letter):
            return True
        grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['d','g','c','e','','a',''],
    ['g','c','','','','',''],
    ['','','f','','','d',''],
    ['e','','','','d','g',''],
    ['','','','d','g','','e'],
    ['a','','','','','','f'],
    ['','','','','e','','a']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# First find diagonal letter if any exists
diag_letter = None
for i in range(7):
    if grid[i][6-i] != '':
        diag_letter = grid[i][6-i]
        break

# Solve
if solve(grid, 0, diag_letter):
    print_grid(grid)
else:
    print("No solution exists")