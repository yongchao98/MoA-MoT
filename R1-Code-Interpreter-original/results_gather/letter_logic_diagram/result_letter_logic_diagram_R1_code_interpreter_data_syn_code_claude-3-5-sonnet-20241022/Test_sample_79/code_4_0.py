def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_empty_cell(grid):
    # First fill minor diagonal cells
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def get_available_letters(grid, row, col):
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(7))
    all_letters = set('abcdefg')
    return all_letters - used_row - used_col

def is_valid(grid, row, col, letter):
    # Minor diagonal must be 'b'
    if row + col == 6 and letter != 'b':
        return False
        
    # Check row
    if letter in grid[row]:
        return False
        
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
        
    return True

def solve(grid):
    empty = get_empty_cell(grid)
    if not empty:
        return True
        
    row, col = empty
    available = get_available_letters(grid, row, col)
    
    # If on minor diagonal, only try 'b'
    if row + col == 6:
        if 'b' not in available:
            return False
        available = {'b'}
    
    for letter in sorted(available):  # Sort for deterministic behavior
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['d','e','c','a','f','','g'],
    ['','c','','','b','g',''],
    ['c','a','f','','','d','e'],
    ['a','','','g','d','e','c'],
    ['f','','g','d','e','c','a'],
    ['b','g','','','c','','f'],
    ['g','d','e','c','','f','']
]

# Verify initial diagonal cells
diag_cells = [(i, 6-i) for i in range(7)]
for row, col in diag_cells:
    if grid[row][col] != '' and grid[row][col] != 'b':
        grid[row][col] = ''  # Clear non-'b' diagonal cells

if solve(grid):
    print_grid(grid)
else:
    print("No solution exists")