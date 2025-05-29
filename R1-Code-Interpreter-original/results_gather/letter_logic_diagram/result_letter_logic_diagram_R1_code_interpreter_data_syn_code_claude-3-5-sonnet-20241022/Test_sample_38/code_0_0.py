def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal if this cell is part of it
    if row + col == 6:  # Cell is on minor diagonal
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
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
    
    # If on minor diagonal, use the first non-empty value on minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
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
grid = [
    ['a','e','','f','g','c',''],
    ['','b','','','','d','a'],
    ['b','f','g','c','d','a','e'],
    ['f','','c','','a','','b'],
    ['g','c','d','a','e','b','f'],
    ['c','','','e','b','','g'],
    ['d','','','','f','g','c']
]

if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")