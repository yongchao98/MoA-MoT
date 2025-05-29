def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_cell_candidates(grid, r, c, initial):
    if initial[r][c] != '':
        return {initial[r][c]}
    
    # If on minor diagonal, must be 'e'
    if r + c == 6:
        return {'e'}
    
    # Get all used values in row and column
    used = set()
    for i in range(7):
        if grid[r][i] != '':
            used.add(grid[r][i])
        if grid[i][c] != '':
            used.add(grid[i][c])
    
    return set('abcdefg') - used

def get_next_cell(grid, initial):
    min_candidates = float('inf')
    best_cell = None
    best_candidates = None
    
    # First check minor diagonal cells
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j, {'e'})
    
    # Then check other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                candidates = get_cell_candidates(grid, i, j, initial)
                if len(candidates) < min_candidates:
                    min_candidates = len(candidates)
                    best_cell = (i, j)
                    best_candidates = candidates
                    if min_candidates == 1:  # Can't get better than this
                        return (i, j, candidates)
    
    return (best_cell[0], best_cell[1], best_candidates) if best_cell else None

def verify_grid(grid):
    # Check rows and columns
    for i in range(7):
        row_vals = [x for x in grid[i] if x != '']
        col_vals = [grid[j][i] for j in range(7) if grid[j][i] != '']
        if len(row_vals) != len(set(row_vals)) or len(col_vals) != len(set(col_vals)):
            return False
    
    # Check minor diagonal
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != 'e':
            return False
    
    return True

def solve(grid, initial):
    if not verify_grid(grid):
        return False
    
    next_cell = get_next_cell(grid, initial)
    if not next_cell:
        return True
    
    r, c, candidates = next_cell
    
    for val in candidates:
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

# Create working grid
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