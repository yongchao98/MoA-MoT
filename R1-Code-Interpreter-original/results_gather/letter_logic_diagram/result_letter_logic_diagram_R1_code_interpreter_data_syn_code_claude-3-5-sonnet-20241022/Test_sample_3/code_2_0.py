def is_valid(grid, row, col, letter, initial_grid):
    # First check: if this is a pre-filled cell, only allow the initial value
    if initial_grid[row][col] != '':
        return letter == initial_grid[row][col]
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def print_grid(grid):
    for row in grid:
        print(','.join(row))

def solve_grid(grid, initial_grid):
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
    letters = 'abcdefg'
    
    # If this is a pre-filled cell, only try that value
    if initial_grid[row][col] != '':
        letters = initial_grid[row][col]
    # If cell is on minor diagonal, match existing diagonal letters
    elif row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['g','','c','a','d','','e'],
    ['','c','','','','',''],
    ['','','','','','g',''],
    ['','','b','','g','f',''],
    ['d','','','','f','','a'],
    ['b','','g','f','','a',''],
    ['','','','','a','d','']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# First, verify all initial values are valid
valid_initial = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            # Temporarily remove this value to check if it's valid
            temp = grid[i][j]
            grid[i][j] = ''
            if not is_valid(grid, i, j, temp, initial_grid):
                valid_initial = False
                break
            grid[i][j] = temp
    if not valid_initial:
        break

if valid_initial and solve_grid(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")