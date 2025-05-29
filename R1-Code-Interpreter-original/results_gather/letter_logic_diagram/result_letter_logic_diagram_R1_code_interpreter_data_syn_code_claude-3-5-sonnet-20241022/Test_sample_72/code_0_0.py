def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal consistency
    if row + col == 6:  # Cell is on minor diagonal
        # Check against first filled cell in minor diagonal
        for i in range(7):
            if grid[i][6-i] and i != row and grid[i][6-i] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, use the first filled cell's letter
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i]:
                letters = grid[i][6-i]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['e','c','b','','f','d','g'],
    ['c','','','f','','','e'],
    ['b','a','','','g','e','c'],
    ['','f','','','e','c',''],
    ['','d','','e','','',''],
    ['','','e','c','b','','f'],
    ['','','c','','','','']
]

# Solve puzzle
if solve(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")