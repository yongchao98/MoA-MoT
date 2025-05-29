def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_move(grid, row, col, letter, initial_grid):
    # Check if this position was pre-filled
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
    
    # If this is on minor diagonal, check if we already have a letter there
    diagonal_letter = None
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
    
    # Try each possible letter
    letters = [diagonal_letter] if diagonal_letter else 'abcdefg'
    for letter in letters:
        if is_valid_move(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'c', '', '', '', 'd', ''],
    ['c', '', '', 'b', 'd', 'f', 'e'],
    ['a', '', '', '', 'f', 'e', ''],
    ['', 'b', '', '', 'e', '', ''],
    ['', '', 'f', '', '', 'a', ''],
    ['d', '', 'e', '', '', '', ''],
    ['', '', '', 'a', 'g', 'b', '']
]

# Create a working copy
grid = [row[:] for row in initial_grid]

# Verify initial state is valid
valid_initial = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            letter = initial_grid[i][j]
            grid[i][j] = ''  # Temporarily remove
            if not is_valid_move(grid, i, j, letter, initial_grid):
                valid_initial = False
            grid[i][j] = letter  # Put back

if valid_initial and solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")