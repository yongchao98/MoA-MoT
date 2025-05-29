def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal constraint
    if row + col == 6:  # Cell is on minor diagonal
        if letter != 'c':
            return False
            
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, only try 'c'
    if row + col == 6:
        if is_valid(grid, row, col, 'c'):
            grid[row][col] = 'c'
            if solve(grid):
                return True
            grid[row][col] = ''
        return False
    
    # Try each letter
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['e','f','','d','','a','c'],
    ['','b','','','a','','e'],
    ['','d','g','a','','e',''],
    ['d','g','a','c','e','f','b'],
    ['g','a','','e','f','b','d'],
    ['','c','','','','d','g'],
    ['','e','f','b','d','g','a']
]

if solve(grid):
    for row in grid:
        print(','.join(row))