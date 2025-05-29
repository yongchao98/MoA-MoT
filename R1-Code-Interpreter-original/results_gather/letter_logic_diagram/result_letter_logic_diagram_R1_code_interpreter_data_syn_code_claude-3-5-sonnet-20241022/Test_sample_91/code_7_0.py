def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check initial grid constraint
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:  # if on minor diagonal
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def get_next_cell(grid):
    # First fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then fill rest of grid row by row
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                return (i, j)
    return None

def solve(grid):
    # Find next empty cell
    pos = get_next_cell(grid)
    if not pos:
        return True
    
    row, col = pos
    
    # If on minor diagonal
    if row + col == 6:
        # Find existing diagonal letter if any
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        
        # Try diagonal letters
        letters = [diag_letter] if diag_letter else list('bcdefg')  # Prefer 'b' for diagonal
        for letter in letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid):
                    return True
                grid[row][col] = ''
    else:
        # Get diagonal letter
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        
        # Try all letters except diagonal letter
        for letter in 'abcdefg':
            if letter != diag_letter and is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid):
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

# First, ensure we can fill the minor diagonal
diag_letter = None
for i in range(7):
    if grid[i][6-i] != '':
        diag_letter = grid[i][6-i]
        break

if not diag_letter:
    # Try 'b' first for diagonal (or any other letter that works)
    diag_letter = 'b'
    # Verify it's possible to use this letter on diagonal
    can_use = True
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != diag_letter:
            can_use = False
            break
    if can_use:
        for i in range(7):
            if grid[i][6-i] == '':
                grid[i][6-i] = diag_letter

# Solve the rest of the grid
if solve(grid):
    print_grid(grid)
else:
    print("No solution exists")