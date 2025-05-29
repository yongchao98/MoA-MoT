def is_valid(grid, row, col, letter, initial_grid):
    # Check if this position was pre-filled
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
    
    # Check minor diagonal - all must be 'g'
    if row + col == 6 and letter != 'g':
        return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, only try 'g'
    if row + col == 6:
        if is_valid(grid, row, col, 'g', initial_grid):
            grid[row][col] = 'g'
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    else:
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial_grid):
                grid[row][col] = letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize the grid with given values
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

# First fill all minor diagonal positions with 'g'
for i in range(7):
    if grid[i][6-i] == '' or grid[i][6-i] == 'g':
        grid[i][6-i] = 'g'

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")