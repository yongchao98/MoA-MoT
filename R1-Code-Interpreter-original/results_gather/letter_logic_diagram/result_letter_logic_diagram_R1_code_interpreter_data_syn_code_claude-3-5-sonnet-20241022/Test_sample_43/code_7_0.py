def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_row_candidates(grid, row, used_in_col):
    # Get all letters that can still be used in this row
    used = set(x for x in grid[row] if x != '')
    candidates = set('abcdefg') - used
    
    # Remove letters that would violate column constraints
    for col in range(7):
        if grid[row][col] == '':
            col_used = used_in_col[col]
            if row + col == 6:  # minor diagonal position
                candidates = {'e'} if 'e' not in used else set()
            else:
                candidates = candidates - col_used
    
    return candidates

def solve_row(grid, row, initial, used_in_col):
    if row == 7:
        return True
    
    # If row is complete, move to next row
    if '' not in grid[row]:
        return solve_row(grid, row + 1, initial, used_in_col)
    
    # Get empty positions in current row
    empty_positions = [(row, j) for j in range(7) if grid[row][j] == '']
    
    # Get all possible candidates for this row
    candidates = get_row_candidates(grid, row, used_in_col)
    if not candidates:
        return False
    
    # Try each empty position
    for r, c in empty_positions:
        # If position is in initial grid, skip it
        if initial[r][c] != '':
            continue
            
        # If on minor diagonal, must be 'e'
        if r + c == 6:
            if 'e' in candidates and 'e' not in used_in_col[c]:
                grid[r][c] = 'e'
                used_in_col[c].add('e')
                if solve_row(grid, row, initial, used_in_col):
                    return True
                grid[r][c] = ''
                used_in_col[c].remove('e')
            continue
            
        # Try each candidate
        for val in candidates:
            if val not in used_in_col[c]:
                grid[r][c] = val
                used_in_col[c].add(val)
                if solve_row(grid, row, initial, used_in_col):
                    return True
                grid[r][c] = ''
                used_in_col[c].remove(val)
                
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

# Create working grid and track used letters in columns
grid = [['' for _ in range(7)] for _ in range(7)]
used_in_col = [set() for _ in range(7)]

# Copy initial values and update used_in_col
for i in range(7):
    for j in range(7):
        if initial[i][j] != '':
            grid[i][j] = initial[i][j]
            used_in_col[j].add(initial[i][j])

# First, ensure all minor diagonal positions are 'e'
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = 'e'
        used_in_col[j].add('e')

# Solve row by row
if solve_row(grid, 0, initial, used_in_col):
    print_grid(grid)
else:
    print("No solution found")