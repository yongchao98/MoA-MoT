def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, prefilled):
    # Must match prefilled value if exists
    if prefilled[row][col] != '' and prefilled[row][col] != letter:
        return False
    
    # Must be 'b' on minor diagonal
    if row + col == 6 and letter != 'b':
        return False
    
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def solve(grid, prefilled, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
        
    # Skip prefilled cells
    if prefilled[row][col] != '':
        grid[row][col] = prefilled[row][col]
        return solve(grid, prefilled, row, col + 1)
    
    # Try each letter
    for letter in ('b' if row + col == 6 else 'abcdefg'):
        if is_valid(grid, row, col, letter, prefilled):
            grid[row][col] = letter
            if solve(grid, prefilled, row, col + 1):
                return True
            grid[row][col] = ''
            
    return False

# Initialize empty grid
grid = [[''] * 7 for _ in range(7)]

# Initialize prefilled values
prefilled = [
    ['c','','a','e','','b',''],
    ['','a','','f','','',''],
    ['','','f','b','g','',''],
    ['','','','','','',''],
    ['','','g','','','',''],
    ['b','','','d','a','e','f'],
    ['','','d','a','','f','b']
]

if solve(grid, prefilled):
    print_grid(grid)
else:
    print("No solution found")