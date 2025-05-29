def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_used_in_row(grid, row):
    return set(x for x in grid[row] if x)

def get_used_in_col(grid, col):
    return set(grid[i][col] for i in range(7) if grid[i][col])

def is_valid_placement(grid, row, col, letter, initial_grid):
    # Check if this is a pre-filled cell
    if initial_grid[row][col] and initial_grid[row][col] != letter:
        return False
        
    # Check if letter already exists in row
    if letter in get_used_in_row(grid, row):
        return False
        
    # Check if letter already exists in column
    if letter in get_used_in_col(grid, col):
        return False
        
    # Check minor diagonal constraint
    if row + col == 6:  # If on minor diagonal
        if initial_grid[row][col]:  # If pre-filled, respect it
            return letter == initial_grid[row][col]
        return letter == 'e'  # Otherwise must be 'e'
            
    return True

def find_next_empty(grid):
    # First find empty cells that have the most constraints
    # (cells with more filled neighbors)
    best_score = -1
    best_pos = None
    
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                score = len(get_used_in_row(grid, i)) + len(get_used_in_col(grid, j))
                if score > best_score:
                    best_score = score
                    best_pos = (i, j)
    
    return best_pos

def solve(grid, initial_grid):
    empty = find_next_empty(grid)
    if not empty:
        return True
        
    row, col = empty
    
    # Get available letters
    used_row = get_used_in_row(grid, row)
    used_col = get_used_in_col(grid, col)
    available = [l for l in 'abcdefg' if l not in used_row and l not in used_col]
    
    # If on minor diagonal and not pre-filled, must be 'e'
    if row + col == 6 and not initial_grid[row][col]:
        available = ['e'] if 'e' in available else []
    
    for letter in available:
        if is_valid_placement(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
            
    return False

def verify_solution(grid):
    # Verify rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Verify minor diagonal
    diag_set = set(grid[i][6-i] for i in range(7))
    if len(diag_set) != 1 or 'e' not in diag_set:
        return False
        
    return True

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
    if verify_solution(grid):
        print_grid(grid)
    else:
        print("Invalid solution found")
else:
    print("No solution exists")