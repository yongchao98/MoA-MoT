def print_grid(grid):
    for row in grid:
        print(','.join(str(x) for x in row))

def copy_grid(grid):
    return [row[:] for row in grid]

def is_valid(grid, row, col, num, initial_grid):
    # Must match initial grid if specified
    if initial_grid[row][col] != '' and initial_grid[row][col] != num:
        return False

    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False

    # Check column    
    for x in range(7):
        if grid[x][col] == num:
            return False

    # Check minor diagonal if this cell is on it
    if row + col == 6:
        diag_val = None
        # Get the first non-empty value on minor diagonal from initial grid
        for i in range(7):
            if initial_grid[i][6-i] != '':
                diag_val = initial_grid[i][6-i]
                break
        if diag_val and num != diag_val:
            return False
        # Check current grid's diagonal
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != num:
                return False

    return True

def find_empty(grid, filled_diagonal):
    # First fill diagonal if not filled
    if not filled_diagonal:
        for i in range(7):
            j = 6 - i
            if grid[i][j] == '':
                return (i, j, False)
    
    # Then fill rest of grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, True)
    return None, None, True

def solve(grid, initial_grid, filled_diagonal=False):
    row, col, new_filled_diagonal = find_empty(grid, filled_diagonal)
    
    if row is None:
        return True

    # If on diagonal, get required letter
    diag_val = None
    if row + col == 6:
        for i in range(7):
            if initial_grid[i][6-i] != '':
                diag_val = initial_grid[i][6-i]
                break

    # Try filling the cell
    if diag_val and row + col == 6:
        # Only try diagonal value if on diagonal
        if is_valid(grid, row, col, diag_val, initial_grid):
            grid[row][col] = diag_val
            if solve(grid, initial_grid, new_filled_diagonal):
                return True
            grid[row][col] = ''
    else:
        # Try all possible values
        for num in 'abcdefg':
            if is_valid(grid, row, col, num, initial_grid):
                grid[row][col] = num
                if solve(grid, initial_grid, new_filled_diagonal):
                    return True
                grid[row][col] = ''

    return False

# Initialize the grid with initial values
initial_grid = [
    ['g','b','e','','d','f',''],
    ['b','','c','','','a',''],
    ['e','c','','f','','','b'],
    ['','d','f','','','',''],
    ['','','','','b','e','c'],
    ['','','','','','','d'],
    ['','','b','','','','']
]

# Create working grid
grid = copy_grid(initial_grid)

# First verify initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            val = initial_grid[i][j]
            initial_grid[i][j] = ''
            if not is_valid(initial_grid, i, j, val, initial_grid):
                valid = False
            initial_grid[i][j] = val
            if not valid:
                break
    if not valid:
        break

print('<<<')
if valid and solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No valid solution exists")
print('>>>')