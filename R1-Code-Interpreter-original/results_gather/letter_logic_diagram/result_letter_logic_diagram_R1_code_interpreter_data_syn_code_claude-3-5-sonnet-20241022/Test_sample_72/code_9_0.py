def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_row_candidates(grid, row, used_letters):
    available = set('abcdefg') - used_letters
    # Remove letters already in this row
    for j in range(7):
        if grid[row][j]:
            available.discard(grid[row][j])
    return available

def get_col_candidates(grid, col, used_letters):
    available = set('abcdefg') - used_letters
    # Remove letters already in this column
    for i in range(7):
        if grid[i][col]:
            available.discard(grid[i][col])
    return available

def find_best_empty(grid):
    min_options = float('inf')
    best_pos = None
    
    # First check cells with most constraints
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                row_filled = sum(1 for x in grid[i] if x)
                col_filled = sum(1 for k in range(7) if grid[k][j])
                constraints = row_filled + col_filled
                if constraints > min_options:
                    min_options = constraints
                    best_pos = (i, j)
    
    return best_pos

def is_valid(grid, row, col, letter, initial_grid):
    # Check initial constraints
    if initial_grid[row][col] and initial_grid[row][col] != letter:
        return False
    
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    return True

def solve(grid, initial_grid, minor_diag_letter):
    pos = find_best_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    
    # Get possible letters
    if row + col == 6:  # On minor diagonal
        candidates = [minor_diag_letter]
    else:
        row_used = set(grid[row])
        col_used = set(grid[i][col] for i in range(7))
        candidates = set('abcdefg') - row_used - col_used
    
    for letter in sorted(candidates):  # Sort for deterministic behavior
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, minor_diag_letter):
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

# Find the letter that should be on the minor diagonal
minor_diag_letter = None
for i in range(7):
    if grid[i][6-i]:
        minor_diag_letter = grid[i][6-i]
        break

if minor_diag_letter:
    # Try to solve with the found minor diagonal letter
    if solve(grid, initial_grid, minor_diag_letter):
        # Verify solution
        valid = True
        # Check rows and columns
        for i in range(7):
            row_set = set(grid[i])
            col_set = set(grid[j][i] for j in range(7))
            if len(row_set) != 7 or len(col_set) != 7:
                valid = False
                break
        
        # Check minor diagonal
        for i in range(7):
            if grid[i][6-i] != minor_diag_letter:
                valid = False
                break
        
        if valid:
            print_grid(grid)
        else:
            print("Invalid solution found")
    else:
        # Try with 'e' if the first attempt failed
        grid = [row[:] for row in initial_grid]
        if solve(grid, initial_grid, 'e'):
            print_grid(grid)
        else:
            print("No solution exists")
else:
    print("No minor diagonal letter found in initial grid")