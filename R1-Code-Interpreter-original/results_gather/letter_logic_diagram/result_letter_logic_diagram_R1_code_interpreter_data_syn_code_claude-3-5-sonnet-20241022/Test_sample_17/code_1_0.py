def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this cell is part of it
    if row + col == 6:  # Cell is on minor diagonal
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid, initial):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and initial[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, we need to use the same letter as other filled cells
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break
        if letters == 'abcdefg':  # No letter found on diagonal yet
            letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial):
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

# Create working grid
grid = [row[:] for row in initial_grid]

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")