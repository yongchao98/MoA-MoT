def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve_grid(grid, initial):
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
    
    # If cell is on minor diagonal, we need to match existing diagonal letters
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, initial):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['g','','c','a','d','','e'],
    ['','c','','','','',''],
    ['','','','','','g',''],
    ['','','b','','g','f',''],
    ['d','','','','f','','a'],
    ['b','','g','f','','a',''],
    ['','','','','a','d','']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# Solve the puzzle
if solve_grid(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")