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
        # Find first filled cell in minor diagonal
        for i in range(7):
            if grid[i][6-i]:
                if letter != grid[i][6-i]:
                    return False
                break
        # If no filled cell found, must be 'e' (from given constraints)
        if letter != 'e':
            return False
            
    return True

def find_next_empty(grid):
    # Find cells on minor diagonal first
    for i in range(7):
        if not grid[i][6-i]:
            return i, 6-i
            
    # Then find other empty cells
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid, initial_grid):
    empty = find_next_empty(grid)
    if not empty:
        return True
        
    row, col = empty
    
    # If on minor diagonal, only try 'e'
    if row + col == 6:
        letters = ['e']
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid_placement(grid, row, col, letter, initial_grid):
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

# First fill all minor diagonal cells with 'e'
for i in range(7):
    if grid[i][6-i] and grid[i][6-i] != 'e':
        print("Contradiction in initial grid")
        exit()
    if not grid[i][6-i]:
        grid[i][6-i] = 'e'

# Solve puzzle
if solve(grid, initial_grid):
    # Verify solution
    valid = True
    for i in range(7):
        row_letters = get_used_in_row(grid, i)
        col_letters = get_used_in_col(grid, i)
        if len(row_letters) != 7 or len(col_letters) != 7:
            valid = False
            break
        if grid[i][6-i] != 'e':
            valid = False
            break
    
    if valid:
        print_grid(grid)
    else:
        print("Invalid solution found")
else:
    print("No solution exists")