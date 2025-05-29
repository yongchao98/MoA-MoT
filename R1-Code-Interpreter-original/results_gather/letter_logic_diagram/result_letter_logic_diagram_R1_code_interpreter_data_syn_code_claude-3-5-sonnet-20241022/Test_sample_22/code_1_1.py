def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    return True

def solve(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True

    # Skip filled positions
    if grid[row][col] != '':
        return solve(grid, initial_grid, row, col + 1)

    # If position was pre-filled in initial grid
    if initial_grid[row][col] != '':
        letter = initial_grid[row][col]
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid, row, col + 1):
                return True
            grid[row][col] = ''
        return False

    # If on minor diagonal
    if row + col == 6:
        if is_valid(grid, row, col, 'c'):
            grid[row][col] = 'c'
            if solve(grid, initial_grid, row, col + 1):
                return True
            grid[row][col] = ''
        return False

    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid, row, col + 1):
                return True
            grid[row][col] = ''

    return False

# Initial grid
grid = [
    ['c','','','','','',''],
    ['','','e','','','',''],
    ['d','e','g','','f','c','b'],
    ['e','g','','f','c','','d'],
    ['','a','f','','','',''],
    ['','f','c','b','','e','g'],
    ['f','c','','','e','g','']
]

# Create working copy and fill diagonal with 'c'
working_grid = [row[:] for row in grid]
for i in range(7):
    j = 6 - i
    if working_grid[i][j] == '':
        working_grid[i][j] = 'c'

if solve(working_grid, grid):
    for row in working_grid:
        print(','.join(row))
else:
    print("No solution exists")