def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_safe(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == num:
            return False
    
    # Check minor diagonal requirement
    if row + col == 6:  # if on minor diagonal
        if num != 'e':
            return False
            
    return True

def find_empty(grid, initial):
    # First fill minor diagonal with 'e'
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j, True)  # True indicates diagonal position
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):  # Skip diagonal positions
                return (i, j, False)
    return None

def solve(grid, initial):
    # Find empty location
    find = find_empty(grid, initial)
    if not find:
        return True
    
    row, col, is_diagonal = find
    
    # If position was pre-filled in initial grid
    if initial[row][col] != '':
        if is_safe(grid, row, col, initial[row][col]):
            grid[row][col] = initial[row][col]
            if solve(grid, initial):
                return True
            grid[row][col] = ''
        return False
    
    # For diagonal positions, only try 'e'
    if is_diagonal:
        if is_safe(grid, row, col, 'e'):
            grid[row][col] = 'e'
            if solve(grid, initial):
                return True
            grid[row][col] = ''
        return False
    
    # Try all possibilities for non-diagonal positions
    for num in 'abcdefg':
        if is_safe(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, initial):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial = [
    ['d','a','','','g','',''],
    ['','','','g','e','',''],
    ['f','','g','','c','',''],
    ['','g','','','','a',''],
    ['g','e','','','','f',''],
    ['e','','','','','b',''],
    ['','d','','f','b','','e']
]

# Create working grid and copy initial values
grid = [['' for _ in range(7)] for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial[i][j] != '':
            grid[i][j] = initial[i][j]

# First, fill all diagonal positions with 'e'
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = 'e'

# Then solve the rest
if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution found")