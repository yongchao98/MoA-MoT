def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_initial_constraints(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

def is_valid(grid, r, c, val, initial, size=7):
    # Check initial constraint
    if initial[r][c] != '' and val != initial[r][c]:
        return False
    
    # Check row
    if val in grid[r]:
        return False
    
    # Check column
    if val in [grid[i][c] for i in range(size) if grid[i][c] != '']:
        return False
    
    # Check minor diagonal constraint
    if r + c == size - 1:  # If cell is on minor diagonal
        # Get the expected diagonal value
        diagonal_val = None
        for i in range(size):
            if initial[i][size-1-i] != '':
                diagonal_val = initial[i][size-1-i]
                break
        if diagonal_val is None:
            diagonal_val = 'e'  # If no pre-filled diagonal value, use 'e'
        if val != diagonal_val:
            return False
    
    return True

def solve(grid, initial, r=0, c=0, size=7):
    if c == size:
        r += 1
        c = 0
    if r == size:
        return check_initial_constraints(grid, initial)
    
    # If this is a pre-filled cell, use it
    if initial[r][c] != '':
        grid[r][c] = initial[r][c]
        return solve(grid, initial, r, c+1)
    
    # Try each letter
    for val in 'abcdefg':
        if is_valid(grid, r, c, val, initial):
            grid[r][c] = val
            if solve(grid, initial, r, c+1):
                return True
            grid[r][c] = ''
    
    return False

# Initialize grid with initial values
initial = [
    ['d','a','','','g','',''],
    ['','','','g','e','',''],
    ['f','','g','','c','',''],
    ['','g','','','','a',''],
    ['g','e','','','','f',''],
    ['e','','','','','b',''],
    ['','d','','f','b','','e']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

# Solve
if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution found")