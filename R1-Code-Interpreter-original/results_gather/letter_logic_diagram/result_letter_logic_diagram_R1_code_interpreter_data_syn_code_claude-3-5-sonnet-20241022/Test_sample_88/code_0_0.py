def is_valid_row(grid, row):
    filled = [x for x in grid[row] if x]
    return len(filled) == len(set(filled))

def is_valid_col(grid, col):
    filled = [grid[r][col] for r in range(7) if grid[r][col]]
    return len(filled) == len(set(filled))

def is_valid_minor_diag(grid):
    diag = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
    return len(set(diag)) <= 1

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def is_valid(grid, row, col, num):
    # Check if number can be placed in current position
    grid[row][col] = num
    valid = (is_valid_row(grid, row) and 
            is_valid_col(grid, col) and 
            is_valid_minor_diag(grid))
    grid[row][col] = ''
    return valid

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, must match other diagonal values
    if row + col == 6:
        diag_vals = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
        if diag_vals:
            letters = diag_vals[0]
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['f','','','','e','a','g'],
    ['c','d','b','','a','g','f'],
    ['','b','','a','g','f','c'],
    ['b','e','','g','f','c','d'],
    ['','a','g','f','c','d','b'],
    ['','g','f','','','b','e'],
    ['g','','','','','e','a']
]

if solve(initial_grid):
    print("<<<")
    for row in initial_grid:
        print(','.join(row))
    print(">>>")
else:
    print("No solution exists")