def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, initial_grid):
    # Must match initial grid if cell was pre-filled
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
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, find the letter that should be there
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if initial_grid[i][j] != '':
                if is_valid(grid, row, col, initial_grid[i][j], initial_grid):
                    grid[row][col] = initial_grid[i][j]
                    if solve_grid(grid, initial_grid):
                        return True
                    grid[row][col] = ''
                return False
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grids
initial_grid = [
    ['g','b','e','','d','f',''],
    ['b','','c','','','a',''],
    ['e','c','','f','','','b'],
    ['','d','f','','','',''],
    ['','','','','b','e','c'],
    ['','','','','','','d'],
    ['','','b','','','','']
]

# Create working grid
grid = [row[:] for row in initial_grid]

print('<<<')
if solve_grid(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")
print('>>>')