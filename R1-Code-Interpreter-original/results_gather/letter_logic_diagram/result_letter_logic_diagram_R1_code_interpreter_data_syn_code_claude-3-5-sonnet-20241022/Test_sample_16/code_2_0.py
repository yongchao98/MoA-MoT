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

def solve(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, get the diagonal letter if exists
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
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
    ['', 'g', '', '', 'c', 'e', 'a'],
    ['g', '', 'b', '', '', '', ''],
    ['', '', 'c', 'e', '', '', ''],
    ['b', '', 'e', 'a', '', 'g', 'd'],
    ['', '', '', '', 'g', '', 'b'],
    ['', '', '', '', 'd', 'b', ''],
    ['', '', '', 'd', 'b', '', '']
]

grid = [row[:] for row in initial_grid]
if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")