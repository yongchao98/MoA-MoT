def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_initial_grid():
    return [
        ['f', 'd', 'a', '.', '.', 'g', '.'],
        ['d', '.', '.', 'b', '.', '.', 'f'],
        ['.', 'c', 'b', '.', 'e', '.', '.'],
        ['.', '.', 'g', 'e', '.', '.', '.'],
        ['.', 'g', 'e', 'f', '.', 'a', 'c'],
        ['.', 'e', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'a', '.', '.', '.']
    ]

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
    if r + c == 6:  # if on minor diagonal
        if val != 'a':
            return False
    
    return True

def find_empty(grid):
    # First check minor diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '.':
            return (i, j)
    
    # Then check remaining positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '.':
                return (i, j)
    return None

def solve(grid, initial):
    # Find empty position
    pos = find_empty(grid)
    if not pos:
        return True
    
    r, c = pos
    
    # If this is a position from initial grid, use that value
    if initial[r][c] != '.':
        if is_valid(grid, r, c, initial[r][c]):
            grid[r][c] = initial[r][c]
            if solve(grid, initial):
                return True
            grid[r][c] = '.'
        return False
    
    # If on minor diagonal, only try 'a'
    if r + c == 6:
        if is_valid(grid, r, c, 'a'):
            grid[r][c] = 'a'
            if solve(grid, initial):
                return True
            grid[r][c] = '.'
        return False
    
    # Try each possible value
    for val in 'abcdefg':
        if is_valid(grid, r, c, val):
            grid[r][c] = val
            if solve(grid, initial):
                return True
            grid[r][c] = '.'
    
    return False

# Initialize grid
initial = get_initial_grid()
grid = [['.' for _ in range(7)] for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial[i][j] != '.':
            grid[i][j] = initial[i][j]

# Solve
if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution exists")