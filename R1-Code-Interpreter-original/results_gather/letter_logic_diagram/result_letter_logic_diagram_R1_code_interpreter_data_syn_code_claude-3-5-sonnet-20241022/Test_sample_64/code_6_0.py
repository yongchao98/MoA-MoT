def is_valid(grid, row, col, letter, initial_grid):
    # Must match initial constraints
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

def solve(grid, initial_grid):
    # Find empty cell
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
    
    # If on minor diagonal, try to use existing diagonal letter
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                if is_valid(grid, row, col, grid[i][6-i], initial_grid):
                    grid[row][col] = grid[i][6-i]
                    if solve(grid, initial_grid):
                        return True
                grid[row][col] = ''
                return False
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
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

# Start with empty grid and copy initial values
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

# First verify initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            temp = initial_grid[i][j]
            initial_grid[i][j] = ''
            if not is_valid(initial_grid, i, j, temp, initial_grid):
                valid = False
            initial_grid[i][j] = temp
            if not valid:
                break

if valid and solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")