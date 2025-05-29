def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check initial constraints
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check minor diagonal constraint
    if row + col == 6:  # if on minor diagonal
        minor_diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                minor_diag_letter = grid[i][6-i]
                break
        if minor_diag_letter and letter != minor_diag_letter:
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then fill rest of grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                return (i, j)
    return None

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    letters = 'abcdefg'
    
    # If on minor diagonal, try to use existing diagonal letter first
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break
    
    # Try each possible letter
    for letter in letters:
        if is_valid(grid, row, col, letter):
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

# First, determine the minor diagonal letter by checking constraints
minor_diag_options = set('abcdefg')
for i in range(7):
    if grid[i][6-i] != '':
        minor_diag_options = {grid[i][6-i]}
        break

# Fill all minor diagonal positions with the same letter first
for letter in minor_diag_options:
    success = True
    temp_grid = [row[:] for row in grid]
    
    # Try to fill minor diagonal
    for i in range(7):
        if temp_grid[i][6-i] == '':
            if is_valid(temp_grid, i, 6-i, letter):
                temp_grid[i][6-i] = letter
            else:
                success = False
                break
    
    if success and solve(temp_grid):
        for row in temp_grid:
            print(','.join(row))
        break