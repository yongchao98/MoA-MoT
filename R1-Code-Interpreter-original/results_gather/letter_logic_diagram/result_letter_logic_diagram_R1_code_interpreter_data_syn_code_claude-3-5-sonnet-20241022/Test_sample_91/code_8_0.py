def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_available_letters(grid, row, col, diag_letter):
    used = set()
    # Check row
    for j in range(7):
        if grid[row][j] != '':
            used.add(grid[row][j])
    
    # Check column
    for i in range(7):
        if grid[i][col] != '':
            used.add(grid[i][col])
    
    # If on minor diagonal, must use diag_letter
    if row + col == 6:
        return {diag_letter} if diag_letter not in used else set()
    
    # Return available letters
    return set('abcdefg') - used

def is_valid_partial(grid):
    # Check initial constraints
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    # Check rows (no duplicates)
    for i in range(7):
        row_letters = [x for x in grid[i] if x != '']
        if len(row_letters) != len(set(row_letters)):
            return False
    
    # Check columns (no duplicates)
    for j in range(7):
        col_letters = [grid[i][j] for i in range(7) if grid[i][j] != '']
        if len(col_letters) != len(set(col_letters)):
            return False
    
    # Check minor diagonal (all same)
    diag_letters = set(grid[i][6-i] for i in range(7) if grid[i][6-i] != '')
    if len(diag_letters) > 1:
        return False
    
    return True

def solve(grid, diag_letter=None):
    if not is_valid_partial(grid):
        return False
    
    # Find empty cell
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
    
    # If this is first diagonal cell, try each possible letter
    if row + col == 6 and diag_letter is None:
        for letter in 'abcdefg':
            available = get_available_letters(grid, row, col, letter)
            if letter in available:
                grid[row][col] = letter
                if solve(grid, letter):
                    return True
                grid[row][col] = ''
    else:
        available = get_available_letters(grid, row, col, diag_letter)
        for letter in available:
            grid[row][col] = letter
            if solve(grid, diag_letter):
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
if solve(grid, diag_letter):
    print_grid(grid)
else:
    print("No solution exists")