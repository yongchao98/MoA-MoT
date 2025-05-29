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
    
    # Check minor diagonal - must be same letter as any existing minor diagonal letter
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
            if initial_grid[i][6-i] != '' and initial_grid[i][6-i] != letter:
                return False
    
    return True

def find_empty(grid):
    # First try minor diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def get_minor_diagonal_letter(initial_grid):
    # Find the letter that appears on the minor diagonal in the initial grid
    for i in range(7):
        if initial_grid[i][6-i] != '':
            return initial_grid[i][6-i]
    return None

def solve(grid, initial_grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    letters = 'abcdefg'
    
    # If on minor diagonal, only try the correct diagonal letter
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(initial_grid)
        if diag_letter:
            letters = diag_letter
    
    for letter in letters:
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'e', '', 'g', 'c', 'f', 'a'],
    ['', '', '', '', '', '', 'b'],
    ['', '', '', 'f', '', '', 'e'],
    ['', '', '', 'a', '', 'e', ''],
    ['', '', '', 'b', '', 'd', ''],
    ['f', '', 'b', 'e', '', '', 'c'],
    ['', '', '', 'd', 'g', '', 'f']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# Get the minor diagonal letter (should be 'a' from position (0,6))
diag_letter = get_minor_diagonal_letter(initial_grid)

# Fill all minor diagonal positions with this letter first
if diag_letter:
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diag_letter

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")