def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_used_in_row(grid, row):
    return set(x for x in grid[row] if x)

def get_used_in_col(grid, col):
    return set(grid[i][col] for i in range(7) if grid[i][col])

def get_candidates(grid, row, col, initial_grid):
    if initial_grid[row][col]:
        return {initial_grid[row][col]}
    
    used_row = get_used_in_row(grid, row)
    used_col = get_used_in_col(grid, col)
    candidates = set('abcdefg') - used_row - used_col
    
    # If on minor diagonal, must be 'e' unless pre-filled
    if row + col == 6:
        if 'e' in candidates:
            return {'e'}
        return set()
        
    return candidates

def find_best_cell(grid):
    min_candidates = float('inf')
    best_pos = None
    best_candidates = None
    
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                candidates = get_candidates(grid, i, j, initial_grid)
                if len(candidates) < min_candidates:
                    min_candidates = len(candidates)
                    best_pos = (i, j)
                    best_candidates = candidates
                    if min_candidates == 1:  # Can't get better than this
                        return best_pos, best_candidates
    
    return best_pos, best_candidates

def verify_partial(grid):
    # Verify no duplicate in rows
    for i in range(7):
        row_values = [x for x in grid[i] if x]
        if len(row_values) != len(set(row_values)):
            return False
    
    # Verify no duplicate in columns
    for j in range(7):
        col_values = [grid[i][j] for i in range(7) if grid[i][j]]
        if len(col_values) != len(set(col_values)):
            return False
    
    # Verify minor diagonal
    diag_values = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
    if diag_values and not all(x == 'e' or x == grid[0][6] for x in diag_values):
        return False
    
    return True

def solve(grid, initial_grid):
    if not verify_partial(grid):
        return False
        
    cell_and_candidates = find_best_cell(grid)
    if not cell_and_candidates[0]:  # No empty cells left
        return True
        
    row, col = cell_and_candidates[0]
    candidates = cell_and_candidates[1]
    
    if not candidates:  # No valid candidates
        return False
    
    for letter in sorted(candidates):  # Sort for deterministic behavior
        grid[row][col] = letter
        if solve(grid, initial_grid):
            return True
        grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['e','c','b','','f','d','g'],
    ['c','','','f','','','e'],
    ['b','a','','','g','e','c'],
    ['','f','','','e','c',''],
    ['','d','','e','','',''],
    ['','','e','c','b','','f'],
    ['','','c','','','','']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# Solve puzzle
if solve(grid, initial_grid):
    # Verify final solution
    all_letters = set('abcdefg')
    valid = True
    
    # Check rows and columns
    for i in range(7):
        if set(grid[i]) != all_letters:
            valid = False
            break
        col_set = {grid[j][i] for j in range(7)}
        if col_set != all_letters:
            valid = False
            break
    
    # Check minor diagonal
    if not all(grid[i][6-i] == 'e' for i in range(7)):
        valid = False
    
    if valid:
        print_grid(grid)
    else:
        print("Invalid solution found")
else:
    print("No solution exists")