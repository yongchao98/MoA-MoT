def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_move(grid, row, col, letter, initial_grid):
    # Check if this position was pre-filled
    if initial_grid[row][col] != '':
        return letter == initial_grid[row][col]
    
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
        return letter == 'c'
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If this is a pre-filled position, use that value
    if initial_grid[row][col] != '':
        if is_valid_move(grid, row, col, initial_grid[row][col], initial_grid):
            grid[row][col] = initial_grid[row][col]
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
        return False
    
    # Try each possible letter
    letters = list('abcdefg')
    # If on minor diagonal, only try 'c'
    if row + col == 6:
        letters = ['c']
        
    for letter in letters:
        if is_valid_move(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize the grid with given values
initial_grid = [
    ['', '', '', '', 'd', 'g', 'c'],
    ['', 'a', '', 'd', '', 'c', ''],
    ['a', '', 'd', 'g', '', '', ''],
    ['', '', 'g', 'c', '', '', 'a'],
    ['d', '', 'c', 'e', '', '', 'b'],
    ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
    ['c', 'e', 'f', 'a', '', 'd', 'g']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# Verify initial grid doesn't violate any constraints
valid_initial = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            temp = initial_grid[i][j]
            initial_grid[i][j] = ''
            if not is_valid_move(initial_grid, i, j, temp, initial_grid):
                valid_initial = False
            initial_grid[i][j] = temp
            if not valid_initial:
                break
    if not valid_initial:
        break

if valid_initial and solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution found")