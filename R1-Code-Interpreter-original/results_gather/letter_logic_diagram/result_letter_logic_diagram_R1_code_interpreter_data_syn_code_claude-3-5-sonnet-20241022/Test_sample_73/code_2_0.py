def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, r, c, val, original, size=7):
    # Check if this position was pre-filled
    if original[r][c] != '' and original[r][c] != val:
        return False
    
    # Check row
    if val in grid[r]:
        return False
    
    # Check column
    if val in [grid[i][c] for i in range(size) if grid[i][c] != '']:
        return False
    
    # Check minor diagonal - must all be the same letter
    if r + c == size - 1:
        diagonal_letter = None
        for i in range(size):
            if grid[i][size-1-i] != '':
                diagonal_letter = grid[i][size-1-i]
                break
        if diagonal_letter and val != diagonal_letter:
            return False
    
    return True

def solve(grid, original, r=0, c=0, size=7):
    if c == size:
        r += 1
        c = 0
    if r == size:
        return True
    
    # Skip if cell is already filled
    if grid[r][c] != '':
        return solve(grid, original, r, c+1)
    
    # Try each letter
    letters = 'abcdefg'
    # If this is on minor diagonal, must match existing diagonal letter
    if r + c == size - 1:
        diagonal_letter = None
        for i in range(size):
            if grid[i][size-1-i] != '':
                diagonal_letter = grid[i][size-1-i]
                break
        if diagonal_letter:
            letters = diagonal_letter
    
    for val in letters:
        if is_valid(grid, r, c, val, original):
            grid[r][c] = val
            if solve(grid, original, r, c+1):
                return True
            grid[r][c] = ''
    
    return False

# Initialize grid and keep original for reference
original = [
    ['f','e','g','a','','',''],
    ['e','g','a','','','','f'],
    ['g','','','','b','f','e'],
    ['','d','','','f','e','g'],
    ['d','','b','','e','','a'],
    ['','b','f','e','','a','d'],
    ['','','e','','a','d','c']
]

grid = [row[:] for row in original]  # Make a copy

# Solve the puzzle
if solve(grid, original):
    print_grid(grid)
else:
    print("No solution found")