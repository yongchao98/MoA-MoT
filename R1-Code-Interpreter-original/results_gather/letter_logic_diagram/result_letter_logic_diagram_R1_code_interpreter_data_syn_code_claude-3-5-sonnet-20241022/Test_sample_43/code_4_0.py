def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, r, c, val):
    # Check row
    for j in range(7):
        if grid[r][j] == val:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][c] == val:
            return False
    
    # Check minor diagonal
    if r + c == 6:  # If on minor diagonal
        if val != 'e':
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid, initial):
    # Find empty location
    empty = find_empty(grid)
    if not empty:
        return True
    
    r, c = empty
    
    # If this is a pre-filled position, use the initial value
    if initial[r][c] != '':
        if is_valid(grid, r, c, initial[r][c]):
            grid[r][c] = initial[r][c]
            if solve(grid, initial):
                return True
            grid[r][c] = ''
        return False
    
    # Try each possible value
    letters = 'e' if r + c == 6 else 'abcdefg'
    for val in letters:
        if is_valid(grid, r, c, val):
            grid[r][c] = val
            if solve(grid, initial):
                return True
            grid[r][c] = ''
    
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

# Solve
if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution found")