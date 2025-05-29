def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, r, c, val, original, size=7):
    # Must match pre-filled position
    if original[r][c] != '' and original[r][c] != val:
        return False
    
    # Check row
    if val in grid[r]:
        return False
    
    # Check column
    if val in [grid[i][c] for i in range(size) if grid[i][c] != '']:
        return False
    
    # Check minor diagonal
    if r + c == size - 1:
        for i in range(size):
            if grid[i][size-1-i] != '' and grid[i][size-1-i] != val:
                return False
    
    return True

def find_empty(grid, original):
    # First fill pre-filled cells
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '' and original[r][c] != '':
                return (r, c)
    # Then fill other cells
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                return (r, c)
    return None

def solve(grid, original):
    empty = find_empty(grid, original)
    if not empty:
        return True
    
    r, c = empty
    letters = 'abcdefg'
    
    # If this is on minor diagonal, must match existing diagonal letter
    if r + c == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break
    
    for val in letters:
        if is_valid(grid, r, c, val, original):
            grid[r][c] = val
            if solve(grid, original):
                return True
            grid[r][c] = ''
    
    return False

# Initialize grid
original = [
    ['f','e','g','a','','',''],
    ['e','g','a','','','','f'],
    ['g','','','','b','f','e'],
    ['','d','','','f','e','g'],
    ['d','','b','','e','','a'],
    ['','b','f','e','','a','d'],
    ['','','e','','a','d','c']
]

grid = [row[:] for row in original]

if solve(grid, original):
    print_grid(grid)
else:
    print("No solution exists")