def is_valid_placement(grid, initial_grid, row, col, letter):
    # Check if this position was pre-filled
    if initial_grid[row][col] != '':
        return letter == initial_grid[row][col]
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, try to find existing diagonal value
    if row + col == 6:
        diag_val = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_val = grid[i][6-i]
                break
        if diag_val:
            if is_valid_placement(grid, initial_grid, row, col, diag_val):
                grid[row][col] = diag_val
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False

    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid_placement(grid, initial_grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize the grid with initial values
initial_grid = [
    ['', 'e', '', 'g', '', 'c', 'b'],
    ['', 'd', '', '', '', 'b', 'a'],
    ['', '', '', 'c', '', '', ''],
    ['', '', 'c', 'b', '', 'e', 'd'],
    ['f', '', '', '', 'e', 'd', ''],
    ['', '', '', '', '', '', ''],
    ['b', 'a', '', 'd', 'g', '', '']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# First, find the letter that should be on the minor diagonal
diag_letter = None
for i in range(7):
    if grid[i][6-i] != '':
        diag_letter = grid[i][6-i]
        break

# Fill in known diagonal values if we found one
if diag_letter:
    for i in range(7):
        if i + (6-i) == 6 and grid[i][6-i] == '':
            grid[i][6-i] = diag_letter

if solve(grid, initial_grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")