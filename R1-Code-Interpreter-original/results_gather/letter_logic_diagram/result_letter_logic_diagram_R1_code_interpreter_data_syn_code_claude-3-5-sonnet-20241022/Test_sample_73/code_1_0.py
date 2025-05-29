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

def solve(grid, r=0, c=0, size=7):
    if c == size:
        r += 1
        c = 0
    if r == size:
        return True
    
    # Skip if cell is already filled
    if grid[r][c] != '':
        return solve(grid, r, c+1)
    
    # Try each letter
    letters = 'abcdefg'
    # If this is the first empty cell in minor diagonal, try all letters
    # Otherwise, must match existing diagonal letter
    if r + c == size - 1:
        diagonal_letter = None
        for i in range(size):
            if grid[i][size-1-i] != '':
                diagonal_letter = grid[i][size-1-i]
                break
        if diagonal_letter:
            letters = diagonal_letter
    
    for val in letters:
        if is_valid(grid, r, c, val):
            grid[r][c] = val
            if solve(grid, r, c+1):
                return True
            grid[r][c] = ''
    
    return False

# Initialize grid
grid = [
    ['f','e','g','a','','',''],
    ['e','g','a','','','','f'],
    ['g','','','','b','f','e'],
    ['','d','','','f','e','g'],
    ['d','','b','','e','','a'],
    ['','b','f','e','','a','d'],
    ['','','e','','a','d','c']
]

# Solve the puzzle
if solve(grid):
    print_grid(grid)
else:
    print("No solution found")