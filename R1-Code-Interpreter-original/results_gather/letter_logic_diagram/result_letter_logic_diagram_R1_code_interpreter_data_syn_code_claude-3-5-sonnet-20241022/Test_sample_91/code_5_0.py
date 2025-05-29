def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, diag_letter):
    # Check initial constraints
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
    if row + col == 6 and letter != diag_letter:
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
            if grid[i][j] == '':
                return (i, j)
    return None

def verify_solution(grid):
    # Check rows and columns
    letters = set('abcdefg')
    for i in range(7):
        if set(grid[i]) != letters:
            return False
        if set(grid[j][i] for j in range(7)) != letters:
            return False
    
    # Check minor diagonal
    diag = set(grid[i][6-i] for i in range(7))
    if len(diag) != 1:
        return False
    
    # Check initial constraints
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    return True

def solve(grid, diag_letter):
    pos = find_empty(grid)
    if not pos:
        return verify_solution(grid)
    
    row, col = pos
    letters = 'abcdefg'
    
    # If on minor diagonal, only try diagonal letter
    if row + col == 6:
        letters = diag_letter
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diag_letter):
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

# Find diagonal letter constraints
diag_letters = set('abcdefg')
for i in range(7):
    if initial_grid[i][6-i] != '':
        diag_letters = {initial_grid[i][6-i]}
        break

# Try each possible diagonal letter
for diag_letter in diag_letters:
    grid = [row[:] for row in initial_grid]
    
    # Pre-fill diagonal
    valid = True
    for i in range(7):
        if grid[i][6-i] == '':
            if is_valid(grid, i, 6-i, diag_letter, diag_letter):
                grid[i][6-i] = diag_letter
            else:
                valid = False
                break
    
    if valid and solve(grid, diag_letter):
        print_grid(grid)
        break