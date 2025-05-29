def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, r, c, val, size=7):
    # Check row
    if val in grid[r]:
        return False
    
    # Check column
    if val in [grid[i][c] for i in range(size) if grid[i][c] != '']:
        return False
    
    # Check minor diagonal constraint
    if r + c == size - 1:  # If cell is on minor diagonal
        # Find first non-empty cell in minor diagonal
        first_val = None
        for i in range(size):
            if grid[i][size-1-i] != '':
                first_val = grid[i][size-1-i]
                break
        if first_val and val != first_val:
            return False
    
    return True

def solve(grid, r=0, c=0, size=7):
    if c == size:
        r += 1
        c = 0
    if r == size:
        return True
    
    if grid[r][c] != '':  # Skip pre-filled cells
        return solve(grid, r, c+1)
    
    # Try each letter
    for val in 'abcdefg':
        if is_valid(grid, r, c, val):
            grid[r][c] = val
            if solve(grid, r, c+1):
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

# Solve
if solve(initial):
    print_grid(initial)
else:
    print("No solution found")