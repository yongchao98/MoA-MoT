def is_valid_placement(grid, row, col, letter, initial_grid):
    # If this position was pre-filled in initial grid, only allow that value
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
    
    # Check minor diagonal consistency
    if row + col == 6:  # If on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve_puzzle(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, try to find existing diagonal value first
    if row + col == 6:
        diagonal_value = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_value = grid[i][6-i]
                break
        if diagonal_value:
            if is_valid_placement(grid, row, col, diagonal_value, initial_grid):
                grid[row][col] = diagonal_value
                if solve_puzzle(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False
    
    for letter in 'abcdefg':
        if is_valid_placement(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_puzzle(grid, initial_grid):
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

# Create a working copy of the grid
grid = [row[:] for row in initial_grid]

# Verify initial grid and solve
if solve_puzzle(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")