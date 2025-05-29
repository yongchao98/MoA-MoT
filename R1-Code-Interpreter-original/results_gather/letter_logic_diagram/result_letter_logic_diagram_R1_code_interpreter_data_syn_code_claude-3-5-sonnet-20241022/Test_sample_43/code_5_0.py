def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_initial_match(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

def get_domain(grid, r, c, initial):
    if initial[r][c] != '':
        return {initial[r][c]}
    if r + c == 6:  # minor diagonal
        return {'e'}
    
    # Get all used values in row and column
    used = set()
    for i in range(7):
        if grid[r][i] != '':
            used.add(grid[r][i])
        if grid[i][c] != '':
            used.add(grid[i][c])
    
    return set('abcdefg') - used

def find_most_constrained(grid, initial):
    min_options = float('inf')
    best_pos = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                domain = get_domain(grid, i, j, initial)
                if len(domain) < min_options:
                    min_options = len(domain)
                    best_pos = (i, j)
    
    return best_pos

def is_valid_so_far(grid):
    # Check rows and columns
    for i in range(7):
        row_vals = [x for x in grid[i] if x != '']
        if len(row_vals) != len(set(row_vals)):
            return False
        
        col_vals = [grid[j][i] for j in range(7) if grid[j][i] != '']
        if len(col_vals) != len(set(col_vals)):
            return False
    
    # Check minor diagonal
    diag_vals = [grid[i][6-i] for i in range(7) if grid[i][6-i] != '']
    if any(x != 'e' for x in diag_vals):
        return False
    
    return True

def solve(grid, initial):
    if not is_valid_so_far(grid):
        return False
        
    pos = find_most_constrained(grid, initial)
    if not pos:
        return True
    
    r, c = pos
    domain = get_domain(grid, r, c, initial)
    
    for val in domain:
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