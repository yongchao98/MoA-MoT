def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal - must be same letter throughout
    if row + col == 6:  # Cell is on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return i, j
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If this is first empty cell on minor diagonal, try all letters
    # Otherwise, must match existing minor diagonal letter
    if row + col == 6:
        minor_diag_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                minor_diag_letter = grid[i][j]
                break
        if minor_diag_letter:
            letters = minor_diag_letter
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'd', '', '', 'b', 'e', 'g'],
    ['', '', '', 'b', '', 'g', ''],
    ['f', '', '', 'e', 'g', '', ''],
    ['a', '', '', '', '', 'd', ''],
    ['', '', 'g', '', '', '', 'a'],
    ['', 'g', '', 'd', '', '', ''],
    ['', 'c', 'd', '', 'a', 'b', '']
]

if solve(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")