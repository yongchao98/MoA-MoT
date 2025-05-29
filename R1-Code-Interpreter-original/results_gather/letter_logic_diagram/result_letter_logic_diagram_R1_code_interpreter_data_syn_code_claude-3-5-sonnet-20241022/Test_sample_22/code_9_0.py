def verify_initial_grid(grid):
    # Check if initial grid has any obvious conflicts
    for i in range(7):
        row_letters = set()
        col_letters = set()
        for j in range(7):
            # Check rows
            if grid[i][j] != '':
                if grid[i][j] in row_letters:
                    return False
                row_letters.add(grid[i][j])
            # Check columns
            if grid[j][i] != '':
                if grid[j][i] in col_letters:
                    return False
                col_letters.add(grid[j][i])
            # Check diagonal
            if i + j == 6 and grid[i][j] != '' and grid[i][j] != 'c':
                return False
    return True

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    # Check minor diagonal
    if row + col == 6 and letter != 'c':
        return False

    return True

def find_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial_grid):
    pos = find_empty(grid)
    if not pos:
        return True

    row, col = pos

    # Get available letters for this position
    if row + col == 6:
        letters = ['c']  # Only 'c' allowed on minor diagonal
    else:
        letters = 'abcdefg'

    # If position was pre-filled, only try that letter
    if initial_grid[row][col] != '':
        letters = [initial_grid[row][col]]

    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''

    return False

# Initial grid
initial_grid = [
    ['c','','','','','',''],
    ['','','e','','','',''],
    ['d','e','g','','f','c','b'],
    ['e','g','','f','c','','d'],
    ['','a','f','','','',''],
    ['','f','c','b','','e','g'],
    ['f','c','','','e','g','']
]

# Create working copy
grid = [row[:] for row in initial_grid]

# First verify initial grid
if not verify_initial_grid(initial_grid):
    print("No solution exists")
else:
    # Fill all empty diagonal positions with 'c'
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = 'c'
    
    # Now solve the rest
    if solve(grid, initial_grid):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution exists")