def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check initial grid constraint
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal (must be same letter)
    if row + col == 6:  # if on minor diagonal
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def get_diagonal_letter():
    # Get the letter that must appear on the minor diagonal
    # by checking pre-filled cells on the diagonal
    for i in range(7):
        if initial_grid[i][6-i] != '':
            return initial_grid[i][6-i]
    return None

def solve(grid, diag_letter):
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
    
    # If this is on minor diagonal
    if row + col == 6:
        if is_valid(grid, row, col, diag_letter):
            grid[row][col] = diag_letter
            if solve(grid, diag_letter):
                return True
            grid[row][col] = ''
    else:
        # Try each letter except diagonal letter
        for letter in 'abcdefg':
            if letter != diag_letter and is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, diag_letter):
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

# First, determine the diagonal letter
diag_letter = get_diagonal_letter()
if not diag_letter:
    # Try each possible letter for diagonal
    for letter in 'abcdefg':
        test_grid = [row[:] for row in grid]
        # Fill diagonal with this letter
        valid = True
        for i in range(7):
            if test_grid[i][6-i] == '':
                if is_valid(test_grid, i, 6-i, letter):
                    test_grid[i][6-i] = letter
                else:
                    valid = False
                    break
        
        if valid and solve(test_grid, letter):
            grid = test_grid
            break
else:
    # Use the pre-determined diagonal letter
    solve(grid, diag_letter)

# Print result
print_grid(grid)