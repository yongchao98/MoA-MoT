def is_valid(grid, row, col, letter, initial_grid):
    # First check: if this position was pre-filled
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
    
    # Check minor diagonal if this position is on it
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve(grid, initial_grid):
    # Find empty position
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, find the letter that should be used
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
        if diagonal_letter:
            if is_valid(grid, row, col, diagonal_letter, initial_grid):
                grid[row][col] = diagonal_letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'g', 'a', '', 'd', '', ''],
    ['', 'a', '', '', '', 'c', 'f'],
    ['a', '', 'd', 'b', 'c', '', 'g'],
    ['e', 'd', '', 'c', 'f', 'g', 'a'],
    ['', 'b', 'c', '', '', 'a', ''],
    ['', 'c', '', 'g', 'a', 'e', 'd'],
    ['', 'f', '', 'a', '', 'd', 'b']
]

# Create a working copy
grid = [row[:] for row in initial_grid]

# Verify initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            letter = initial_grid[i][j]
            grid[i][j] = ''
            if not is_valid(grid, i, j, letter, initial_grid):
                valid = False
                break
            grid[i][j] = letter
    if not valid:
        break

if valid and solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")