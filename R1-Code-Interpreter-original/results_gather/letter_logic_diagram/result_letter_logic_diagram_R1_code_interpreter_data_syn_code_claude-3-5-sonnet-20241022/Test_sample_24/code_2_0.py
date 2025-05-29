def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_initial_constraints(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '.' and grid[i][j] != initial[i][j]:
                return False
    return True

def is_valid(grid, r, c, val, initial, size=7):
    # Check if we're violating initial constraints
    if initial[r][c] != '.' and val != initial[r][c]:
        return False
    
    # Check row
    if val in grid[r]:
        return False
    
    # Check column
    if val in [grid[i][c] for i in range(size) if grid[i][c] != '.']:
        return False
    
    # Check minor diagonal
    if r + c == size - 1:
        # Find the first non-empty value on minor diagonal
        diag_val = None
        for i in range(size):
            j = size - 1 - i
            if initial[i][j] != '.':
                diag_val = initial[i][j]
                break
        if diag_val and val != diag_val:
            return False
        
    return True

def solve(grid, initial, r=0, c=0, size=7):
    if c == size:
        r += 1
        c = 0
    if r == size:
        return check_initial_constraints(grid, initial)
    
    # If cell is pre-filled in initial grid, use that value
    if initial[r][c] != '.':
        grid[r][c] = initial[r][c]
        return solve(grid, initial, r, c+1)
    
    # Find diagonal value if we're on diagonal
    diag_val = None
    if r + c == size - 1:
        for i in range(size):
            j = size - 1 - i
            if initial[i][j] != '.':
                diag_val = initial[i][j]
                break
    
    # Try each possible letter
    for val in ['a'] if diag_val == 'a' and r + c == size - 1 else 'abcdefg':
        if is_valid(grid, r, c, val, initial):
            grid[r][c] = val
            if solve(grid, initial, r, c+1):
                return True
            grid[r][c] = '.'
    return False

# Initialize grid with the given values
initial_grid = [
    ['f', 'd', 'a', '.', '.', 'g', '.'],
    ['d', '.', '.', 'b', '.', '.', 'f'],
    ['.', 'c', 'b', '.', 'e', '.', '.'],
    ['.', '.', 'g', 'e', '.', '.', '.'],
    ['.', 'g', 'e', 'f', '.', 'a', 'c'],
    ['.', 'e', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'a', '.', '.', '.']
]

# Create working grid
working_grid = [['.' for _ in range(7)] for _ in range(7)]

if solve(working_grid, initial_grid):
    print_grid(working_grid)
else:
    print("No solution exists")