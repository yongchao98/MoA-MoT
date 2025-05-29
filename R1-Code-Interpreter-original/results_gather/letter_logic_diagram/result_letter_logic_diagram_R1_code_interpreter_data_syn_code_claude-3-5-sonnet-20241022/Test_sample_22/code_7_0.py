def is_valid_move(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    # If on minor diagonal, must be 'c'
    if row + col == 6 and letter != 'c':
        return False

    return True

def find_next_empty(grid):
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

def solve_puzzle(grid, initial_grid):
    # Find next empty position
    pos = find_next_empty(grid)
    if not pos:
        return True

    row, col = pos

    # If this position was in initial grid, must use that value
    if initial_grid[row][col] != '':
        if is_valid_move(grid, row, col, initial_grid[row][col]):
            grid[row][col] = initial_grid[row][col]
            if solve_puzzle(grid, initial_grid):
                return True
            grid[row][col] = ''
        return False

    # If on minor diagonal, must be 'c'
    if row + col == 6:
        if is_valid_move(grid, row, col, 'c'):
            grid[row][col] = 'c'
            if solve_puzzle(grid, initial_grid):
                return True
            grid[row][col] = ''
        return False

    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter):
            grid[row][col] = letter
            if solve_puzzle(grid, initial_grid):
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

# First verify initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            temp = initial_grid[i][j]
            initial_grid[i][j] = ''
            if not is_valid_move(grid, i, j, temp):
                valid = False
            initial_grid[i][j] = temp

# Fill in all diagonal positions with 'c' first
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = 'c'

if solve_puzzle(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")