def is_valid_position(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
            
    # Check minor diagonal
    if row + col == 6:  # if on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
                
    return True

def solve_grid(grid, initial_grid, pos=0):
    if pos == 49:  # 7x7 = 49 positions
        return True
        
    row = pos // 7
    col = pos % 7
    
    # Skip pre-filled positions
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve_grid(grid, initial_grid, pos + 1)
    
    # If on minor diagonal, must be 'g'
    if row + col == 6:
        if is_valid_position(grid, row, col, 'g'):
            grid[row][col] = 'g'
            if solve_grid(grid, initial_grid, pos + 1):
                return True
        grid[row][col] = ''
        return False
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid_position(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''
            
    return False

# Initialize grid with initial positions
initial_grid = [
    ['', 'g', 'f', '', 'e', '', ''],
    ['', 'f', 'a', '', 'c', '', 'b'],
    ['', '', '', 'c', 'd', '', ''],
    ['', '', 'c', 'd', 'b', '', ''],
    ['e', 'c', '', '', '', '', 'a'],
    ['', 'd', '', 'g', 'f', 'a', 'e'],
    ['', '', 'g', 'f', '', 'e', 'c']
]

grid = [['' for _ in range(7)] for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

if solve_grid(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")