def is_valid(grid, row, col, letter, initial_grid):
    # Check if this position was pre-filled
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this position is on it
    if row + col == 6:  # Position is on minor diagonal
        # Check if it matches the first filled cell on minor diagonal
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
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
    letters = 'abcdefg'
    
    # If position is on minor diagonal, we need to match existing diagonal letter
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break
    
    for letter in letters:
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

grid = [row[:] for row in initial_grid]  # Make a copy

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")