def print_grid(grid):
    for row in grid:
        print(','.join(row))

def verify_partial(grid, r, c, val):
    # Check row
    for j in range(7):
        if j != c and grid[r][j] == val:
            return False
    
    # Check column
    for i in range(7):
        if i != r and grid[i][c] == val:
            return False
    
    # Check minor diagonal
    if r + c == 6:
        if val != 'e':
            return False
    
    return True

def get_candidates(grid, r, c, initial):
    if initial[r][c] != '':
        return [initial[r][c]]
    if r + c == 6:
        return ['e']
    
    candidates = []
    for val in 'abcdefg':
        if verify_partial(grid, r, c, val):
            candidates.append(val)
    return candidates

def find_best_cell(grid, initial):
    min_candidates = 8  # More than possible (7)
    best_pos = None
    best_candidates = None
    
    # First priority: minor diagonal cells
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j, ['e'])
    
    # Second priority: cells with fewest candidates
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                candidates = get_candidates(grid, i, j, initial)
                if len(candidates) < min_candidates:
                    min_candidates = len(candidates)
                    best_pos = (i, j)
                    best_candidates = candidates
                    if min_candidates == 1:  # Can't get better than this
                        return (i, j, candidates)
    
    return (best_pos[0], best_pos[1], best_candidates) if best_pos else None

def solve(grid, initial):
    # Find best empty cell
    cell = find_best_cell(grid, initial)
    if not cell:
        return True  # Grid is complete
    
    r, c, candidates = cell
    
    # Try each candidate
    for val in candidates:
        if verify_partial(grid, r, c, val):
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

# Create working grid and copy initial values
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial[i][j] != '':
            if not verify_partial(grid, i, j, initial[i][j]):
                print("Initial configuration is invalid")
                exit()
            grid[i][j] = initial[i][j]

# Pre-fill minor diagonal with 'e' where empty
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = 'e'

if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution found")