def get_diagonal_letter(grid):
    # Check existing letters on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid_placement(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['g','b','f','','c','a',''],
    ['','f','d','','a','e',''],
    ['','d','c','','','',''],
    ['','','','','','',''],
    ['','','e','','b','',''],
    ['a','e','','','','d','c'],
    ['e','','','','','c','']
]

# First, find and fill the minor diagonal
diag_letter = get_diagonal_letter(grid=initial_grid)
if diag_letter is None:
    # If no letter is found on diagonal, try each possibility
    for letter in 'abcdefg':
        test_grid = [row[:] for row in initial_grid]
        valid = True
        # Try to fill diagonal with this letter
        for i in range(7):
            if test_grid[i][6-i] != '' and test_grid[i][6-i] != letter:
                valid = False
                break
            test_grid[i][6-i] = letter
        if valid:
            diag_letter = letter
            initial_grid = test_grid
            break
else:
    # Fill all diagonal positions with the found letter
    for i in range(7):
        initial_grid[i][6-i] = diag_letter

# Now solve the rest of the grid
if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")