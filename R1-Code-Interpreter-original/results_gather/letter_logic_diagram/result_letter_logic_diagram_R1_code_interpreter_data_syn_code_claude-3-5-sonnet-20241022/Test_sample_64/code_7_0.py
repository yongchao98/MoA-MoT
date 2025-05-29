def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal - must be 'a'
    if row + col == 6 and letter != 'a':
        return False
    
    return True

def solve(grid, initial_grid):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (initial_grid[i][j] == '' or initial_grid[i][j] == grid[i][j]):
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    
    # If this position was pre-filled in initial grid, use that value
    if initial_grid[row][col] != '':
        if is_valid(grid, row, col, initial_grid[row][col]):
            grid[row][col] = initial_grid[row][col]
            if solve(grid, initial_grid):
                return True
        grid[row][col] = ''
        return False
    
    # If on minor diagonal, must use 'a'
    if row + col == 6:
        if is_valid(grid, row, col, 'a'):
            grid[row][col] = 'a'
            if solve(grid, initial_grid):
                return True
        grid[row][col] = ''
        return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['g', '', '', 'd', '', '', ''],
    ['', 'b', 'd', '', '', '', ''],
    ['', '', '', '', 'e', 'g', 'a'],
    ['', '', 'f', '', '', 'a', ''],
    ['', 'f', 'e', '', '', 'b', ''],
    ['f', '', '', '', 'b', 'd', 'c'],
    ['', 'g', 'a', '', '', '', 'f']
]

# Create working grid and copy initial values
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

# Fill all minor diagonal positions with 'a'
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = 'a'

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")