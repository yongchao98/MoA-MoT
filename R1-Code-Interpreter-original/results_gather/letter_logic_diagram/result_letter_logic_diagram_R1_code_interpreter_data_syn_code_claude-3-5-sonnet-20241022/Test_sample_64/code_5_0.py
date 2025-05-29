def verify_initial_constraints(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def is_valid(grid, row, col, letter, initial_grid):
    # Must match initial grid if cell was pre-filled
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
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid, initial_grid, pos=0):
    if pos == 49:
        return verify_initial_constraints(grid, initial_grid)
    
    row = pos // 7
    col = pos % 7
    
    # If this is a pre-filled position, use the initial value
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, pos + 1)
    
    # If on minor diagonal, must match existing diagonal letter if any
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter:
            if is_valid(grid, row, col, diagonal_letter, initial_grid):
                grid[row][col] = diagonal_letter
                if solve(grid, initial_grid, pos + 1):
                    return True
            grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, pos + 1):
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

# Start with empty grid
grid = [[''] * 7 for _ in range(7)]

# First copy initial values
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")