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

def is_valid(grid, r, c, val, initial):
    # Check initial constraints
    if initial[r][c] != '.' and initial[r][c] != val:
        return False

    # Check row
    if val in grid[r]:
        return False

    # Check column
    for i in range(7):
        if grid[i][c] == val:
            return False

    # Check minor diagonal requirement
    if r + c == 6:  # if on minor diagonal
        # Must be 'a' on minor diagonal
        if val != 'a':
            return False
    
    return True

def find_empty(grid):
    # Find empty cell, prioritizing minor diagonal
    # First check minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '.':
            return (i, j)
    
    # Then check rest of grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '.':
                return (i, j)
    return None

def solve(grid, initial):
    pos = find_empty(grid)
    if not pos:
        return True
    
    r, c = pos
    letters = ['a'] if r + c == 6 else 'abcdefg'
    
    for val in letters:
        if is_valid(grid, r, c, val, initial):
            # Try placing the value
            grid[r][c] = val
            
            # Recursively try to solve the rest
            if solve(grid, initial):
                return True
            
            # If we get here, we need to backtrack
            grid[r][c] = '.'
    
    return False

# Initialize grid
initial = get_initial_grid()
grid = [row[:] for row in initial]

# First, ensure all minor diagonal positions are 'a'
for i in range(7):
    j = 6 - i
    if grid[i][j] == '.':
        grid[i][j] = 'a'
    elif grid[i][j] != 'a':
        print("Invalid initial grid - minor diagonal must be 'a'")
        exit()

# Now solve the rest
if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution exists")