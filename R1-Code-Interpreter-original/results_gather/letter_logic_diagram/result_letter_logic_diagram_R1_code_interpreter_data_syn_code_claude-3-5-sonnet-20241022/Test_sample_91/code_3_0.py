def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, initial_grid):
    # Check if this conflicts with initial grid
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
    
    # Check minor diagonal (top-right to bottom-left)
    if row + col == 6:  # If this cell is on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def get_minor_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def find_empty(grid):
    # Find next empty cell, prioritizing minor diagonal
    # First check minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    # Then check rest of grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):  # Skip minor diagonal
                return (i, j)
    return None

def solve(grid, initial_grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    
    # If we're on minor diagonal
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter:
            # Must use the same letter as rest of diagonal
            if is_valid(grid, row, col, diag_letter, initial_grid):
                grid[row][col] = diag_letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
        else:
            # First diagonal cell - try each letter
            for letter in 'abcdefg':
                if is_valid(grid, row, col, letter, initial_grid):
                    grid[row][col] = letter
                    if solve(grid, initial_grid):
                        return True
                    grid[row][col] = ''
    else:
        # Regular cell - try each letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial_grid):
                grid[row][col] = letter
                if solve(grid, initial_grid):
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

# Solve
if solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")