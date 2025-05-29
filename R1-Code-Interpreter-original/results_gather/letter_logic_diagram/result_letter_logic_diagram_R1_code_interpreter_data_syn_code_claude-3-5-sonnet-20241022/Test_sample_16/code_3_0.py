def get_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter, initial_grid):
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
        diag_letter = get_diagonal_letter(grid)
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def find_empty(grid):
    # First try to fill minor diagonal
    diag_letter = get_diagonal_letter(grid)
    if diag_letter:
        for i in range(7):
            j = 6 - i
            if grid[i][j] == '':
                return (i, j)
    
    # If no empty cells on diagonal or no diagonal letter yet,
    # find first empty cell on diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Then fill other cells
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
    
    # If on minor diagonal
    if row + col == 6:
        diag_letter = get_diagonal_letter(grid)
        letters = [diag_letter] if diag_letter else 'bcdefg'  # Using 'b' as primary choice
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'g', '', '', 'c', 'e', 'a'],
    ['g', '', 'b', '', '', '', ''],
    ['', '', 'c', 'e', '', '', ''],
    ['b', '', 'e', 'a', '', 'g', 'd'],
    ['', '', '', '', 'g', '', 'b'],
    ['', '', '', '', 'd', 'b', ''],
    ['', '', '', 'd', 'b', '', '']
]

grid = [row[:] for row in initial_grid]

# First ensure we have 'b' on the minor diagonal where pre-filled
for i in range(7):
    j = 6 - i
    if initial_grid[i][j] != '' and initial_grid[i][j] != 'b':
        print("No solution exists - diagonal constraint cannot be satisfied")
        exit()

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")