def print_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def check_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] and grid[i][j] != initial[i][j]:
                return False
    return True

def is_valid(grid, row, col, letter, initial):
    # Check initial constraint
    if initial[row][col] and initial[row][col] != letter:
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
            j = 6 - i
            if grid[i][j] and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid, initial):
    # First, fill cells that are in initial grid
    for i in range(7):
        for j in range(7):
            if initial[i][j] and not grid[i][j]:
                return (i, j)
    
    # Then, fill minor diagonal cells
    for i in range(7):
        j = 6 - i
        if not grid[i][j]:
            return (i, j)
    
    # Finally, fill remaining cells
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return (i, j)
    
    return None

def solve(grid, initial):
    pos = find_empty(grid, initial)
    if not pos:
        return True
    
    row, col = pos
    
    # If this is an initial position, only try that letter
    if initial[row][col]:
        letter = initial[row][col]
        if is_valid(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial):
                return True
            grid[row][col] = ''
        return False
    
    # If on minor diagonal, try to match existing diagonal letter
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if grid[i][6-i]:
                diag_letter = grid[i][6-i]
                break
        if diag_letter:
            if is_valid(grid, row, col, diag_letter, initial):
                grid[row][col] = diag_letter
                if solve(grid, initial):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial):
                return True
            grid[row][col] = ''
    
    return False

# Initialize puzzle
initial = [
    ['c','','g','','f','e',''],
    ['a','g','','','','d',''],
    ['','','','','d','',''],
    ['','','','d','','a','g'],
    ['','e','','c','a','g','b'],
    ['e','','','a','g','b','f'],
    ['d','','a','','b','','']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

# First pass: fill in initial values
for i in range(7):
    for j in range(7):
        if initial[i][j]:
            grid[i][j] = initial[i][j]

if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution exists")