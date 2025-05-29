def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
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
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

# Try each possible letter for the minor diagonal
for diagonal_letter in 'abcdefg':
    grid = [row[:] for row in initial_grid]  # Make a copy
    
    # Check if this diagonal letter is compatible with existing positions
    valid = True
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '' and grid[i][j] != diagonal_letter:
            valid = False
            break
    
    if not valid:
        continue
        
    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = diagonal_letter
    
    # Try to solve the rest of the grid
    if solve_grid(grid):
        for row in grid:
            print(','.join(row))
        break