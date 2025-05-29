def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_complete(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if not all(x == diag[0] for x in diag):
        return False
    
    return True

def is_valid_partial(grid, r, c, val, initial):
    # Check if we're violating initial constraints
    if initial[r][c] != '.' and val != initial[r][c]:
        return False
    
    # Check row
    row = [x for x in grid[r] if x != '.']
    if val in row:
        return False
    
    # Check column
    col = [grid[i][c] for i in range(7) if grid[i][c] != '.']
    if val in col:
        return False
    
    # Check minor diagonal
    if r + c == 6:  # if on minor diagonal
        diag_vals = []
        for i in range(7):
            if grid[i][6-i] != '.':
                diag_vals.append(grid[i][6-i])
        if diag_vals and val != diag_vals[0]:
            return False
    
    return True

def solve(grid, initial, r=0, c=0):
    if c == 7:
        r += 1
        c = 0
    if r == 7:
        return is_valid_complete(grid)
    
    # If cell is pre-filled in initial grid
    if initial[r][c] != '.':
        grid[r][c] = initial[r][c]
        return solve(grid, initial, r, c+1)
    
    # If on minor diagonal, find the required value
    diag_val = None
    if r + c == 6:
        for i in range(7):
            if grid[i][6-i] != '.':
                diag_val = grid[i][6-i]
                break
        if diag_val is None:
            for i in range(7):
                if initial[i][6-i] != '.':
                    diag_val = initial[i][6-i]
                    break
    
    # Try values
    candidates = [diag_val] if diag_val else 'abcdefg'
    for val in candidates:
        if is_valid_partial(grid, r, c, val, initial):
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