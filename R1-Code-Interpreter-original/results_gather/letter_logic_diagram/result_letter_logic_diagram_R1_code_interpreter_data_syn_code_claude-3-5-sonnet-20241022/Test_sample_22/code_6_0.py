def is_valid(grid, row, col, letter, initial_grid):
    # Must respect initial grid values
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    # Check minor diagonal - must be 'c'
    if row + col == 6 and letter != 'c':
        return False

    return True

def solve(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True

    # If this position is pre-filled in initial grid, move to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, row, col + 1)

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, row, col + 1):
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

# First verify that initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            temp = initial_grid[i][j]
            initial_grid[i][j] = ''
            if not is_valid(grid, i, j, temp, initial_grid):
                valid = False
            initial_grid[i][j] = temp
            if i + j == 6 and temp != 'c':  # verify diagonal is 'c'
                valid = False

if valid and solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")