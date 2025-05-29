def verify_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] and grid[i][j] != initial[i][j]:
                return False
    return True

def verify_solution(grid):
    # Check rows
    for row in grid:
        if set(row) != set('abcdefg'):
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if set(col) != set('abcdefg'):
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:
        return False
    
    return True

def get_next_cell(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def get_valid_letters(grid, row, col):
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(7))
    
    # If on minor diagonal, must match existing diagonal letter
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i]:
                return {grid[i][6-i]}
    
    return set('abcdefg') - used_row - used_col

def solve(grid, initial):
    # Find empty cell
    next_cell = get_next_cell(grid)
    if not next_cell:
        return verify_solution(grid)
    
    row, col = next_cell
    
    # If this is a pre-filled cell, use the initial value
    if initial[row][col]:
        grid[row][col] = initial[row][col]
        if solve(grid, initial):
            return True
        grid[row][col] = ''
        return False
    
    valid_letters = get_valid_letters(grid, row, col)
    for letter in valid_letters:
        grid[row][col] = letter
        if solve(grid, initial):
            return True
        grid[row][col] = ''
    
    return False

# Initialize the puzzle
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

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial[i][j]:
            grid[i][j] = initial[i][j]

if solve(grid, initial):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")