def verify_initial_values(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def is_valid_placement(grid, initial_grid, row, col, letter):
    # Must match initial value if one exists
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

def find_next_cell(grid):
    # First try to fill minor diagonal
    for i in range(7):
        if i + (6-i) == 6 and grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial_grid):
    # Find next empty cell
    next_cell = find_next_cell(grid)
    if not next_cell:
        return verify_initial_values(grid, initial_grid)
    
    row, col = next_cell
    
    # If on minor diagonal, use the diagonal letter
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        
        if diag_letter:
            if is_valid_placement(grid, initial_grid, row, col, diag_letter):
                grid[row][col] = diag_letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid_placement(grid, initial_grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize the grid
initial_grid = [
    ['', 'e', '', 'g', '', 'c', 'b'],
    ['', 'd', '', '', '', 'b', 'a'],
    ['', '', '', 'c', '', '', ''],
    ['', '', 'c', 'b', '', 'e', 'd'],
    ['f', '', '', '', 'e', 'd', ''],
    ['', '', '', '', '', '', ''],
    ['b', 'a', '', 'd', 'g', '', '']
]

# Create working grid and copy initial values
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

# Find diagonal letter
diag_letter = None
for i in range(7):
    if initial_grid[i][6-i] != '':
        diag_letter = initial_grid[i][6-i]
        break

# Pre-fill diagonal if we know the letter
if diag_letter:
    for i in range(7):
        if i + (6-i) == 6 and grid[i][6-i] == '':
            grid[i][6-i] = diag_letter

print('<<<')
if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')