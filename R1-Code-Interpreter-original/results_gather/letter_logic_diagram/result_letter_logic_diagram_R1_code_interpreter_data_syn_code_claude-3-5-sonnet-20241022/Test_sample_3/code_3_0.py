def is_valid_placement(grid, row, col, letter, initial_grid):
    # If this position has an initial value, it must match
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
            
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
                
    return True

def initialize_grid():
    return [
        ['g','','c','a','d','','e'],
        ['','c','','','','',''],
        ['','','','','','g',''],
        ['','','b','','g','f',''],
        ['d','','','','f','','a'],
        ['b','','g','f','','a',''],
        ['','','','','a','d','']
    ]

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
    
    # If on minor diagonal, find the letter that should be there
    diagonal_letter = None
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
    
    # Try each possible letter
    for letter in 'abcdefg':
        # If we're on minor diagonal and have a letter, only try that
        if row + col == 6 and diagonal_letter and letter != diagonal_letter:
            continue
            
        if is_valid_placement(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            
            if solve(grid, initial_grid):
                return True
                
            grid[row][col] = ''
            
    return False

# Initialize grids
initial_grid = initialize_grid()
working_grid = [row[:] for row in initial_grid]

# First verify all initial positions are valid
valid_initial = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            temp = working_grid[i][j]
            working_grid[i][j] = ''
            if not is_valid_placement(working_grid, i, j, temp, initial_grid):
                valid_initial = False
                break
            working_grid[i][j] = temp

if not valid_initial:
    print("Initial configuration is invalid")
else:
    # Find the letter that should be on the minor diagonal
    diagonal_letter = None
    for i in range(7):
        if working_grid[i][6-i] != '':
            diagonal_letter = working_grid[i][6-i]
            break
    
    # Fill in the minor diagonal first if we know what it should be
    if diagonal_letter:
        for i in range(7):
            if working_grid[i][6-i] == '':
                working_grid[i][6-i] = diagonal_letter

    if solve(working_grid, initial_grid):
        for row in working_grid:
            print(','.join(row))
    else:
        print("No solution exists")