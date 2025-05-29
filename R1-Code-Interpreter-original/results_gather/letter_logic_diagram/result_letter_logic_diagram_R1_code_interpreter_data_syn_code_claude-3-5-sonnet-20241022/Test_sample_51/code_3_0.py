def is_valid(grid, row, col, letter, initial_grid):
    # Check if this position was pre-filled
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
    
    # Check minor diagonal constraint
    if row + col == 6:  # Cell is on minor diagonal
        if letter != 'd':  # Must be 'd' on minor diagonal
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        j = 6-i
        if grid[i][j] == '':
            return i, j
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = list('abcdefg')
    
    # If cell is on minor diagonal, only try 'd'
    if row + col == 6:
        letters = ['d']
    
    # Try each possible letter
    for letter in letters:
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['b','','','','g','e','d'],
    ['','','','','','d','b'],
    ['','','','e','','',''],
    ['','','','','b','','f'],
    ['','','d','','','f',''],
    ['e','d','','','f','',''],
    ['d','b','','f','','g','e']
]

# Create a deep copy to work with
grid = [row[:] for row in initial_grid]

# First, fill all minor diagonal cells with 'd'
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = 'd'

# Solve the puzzle
if solve_grid(grid, initial_grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")