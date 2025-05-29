def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:  # If this is on minor diagonal
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    # Check if this placement matches initial grid
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    return True

def solve(grid, minor_diag_letter=None):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    
    # If this is the first cell on minor diagonal we're filling
    if row + col == 6 and minor_diag_letter is None:
        # Try each letter for the entire diagonal
        for letter in 'abcdefg':
            if all(grid[i][6-i] == '' or grid[i][6-i] == letter 
                  for i in range(7)):
                minor_diag_letter = letter
                if try_letter(grid, row, col, letter, minor_diag_letter):
                    return True
    else:
        # If on minor diagonal, must use minor_diag_letter
        if row + col == 6:
            if try_letter(grid, row, col, minor_diag_letter, minor_diag_letter):
                return True
        else:
            # Try all letters
            for letter in 'abcdefg':
                if try_letter(grid, row, col, letter, minor_diag_letter):
                    return True
    
    return False

def try_letter(grid, row, col, letter, minor_diag_letter):
    if is_valid(grid, row, col, letter):
        grid[row][col] = letter
        if solve(grid, minor_diag_letter):
            return True
        grid[row][col] = ''
    return False

# Initial grid
initial_grid = [
    ['d','g','c','e','','a',''],
    ['g','c','','','','',''],
    ['','','f','','','d',''],
    ['e','','','','d','g',''],
    ['','','','d','g','','e'],
    ['a','','','','','','f'],
    ['','','','','e','','a']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# Solve and print
if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")