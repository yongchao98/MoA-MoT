def is_valid(grid, row, col, num):
    # Check if number is valid in row
    if num in grid[row]:
        return False
    
    # Check if number is valid in column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    # If on diagonal, must be 'e'
    if row + col == 6 and num != 'e':
        return False
        
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial):
    find = find_empty(grid)
    if not find:
        return True
    
    row, col = find
    
    # If this is a pre-filled position, skip it
    if initial[row][col] != '':
        grid[row][col] = initial[row][col]
        return solve(grid, initial)
    
    # If this is a diagonal position, must be 'e'
    if row + col == 6:
        if is_valid(grid, row, col, 'e'):
            grid[row][col] = 'e'
            if solve(grid, initial):
                return True
        grid[row][col] = ''
        return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial):
                return True
            grid[row][col] = ''
    
    return False

# Initialize empty grid
grid = [[''] * 7 for _ in range(7)]

# Initial configuration
initial = [
    ['g','b','f','','c','a',''],
    ['','f','d','','a','e',''],
    ['','d','c','','','',''],
    ['','','','','','',''],
    ['','','e','','b','',''],
    ['a','e','','','','d','c'],
    ['e','','','','','c','']
]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial[i][j] != '':
            grid[i][j] = initial[i][j]

# Fill diagonal with 'e' first
for i in range(7):
    grid[i][6-i] = 'e'

# Solve the puzzle
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")