def is_valid_row(grid, row):
    filled = [x for x in grid[row] if x]
    return len(filled) == len(set(filled))

def is_valid_col(grid, col):
    filled = [grid[r][col] for r in range(7) if grid[r][col]]
    return len(filled) == len(set(filled))

def is_valid_minor_diagonal(grid):
    diag = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
    return len(set(diag)) <= 1

def matches_initial(grid, initial):
    return all(grid[i][j] == initial[i][j] for i in range(7) for j in range(7) if initial[i][j])

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def get_valid_letters(grid, row, col, initial):
    if initial[row][col]:
        return [initial[row][col]]
    
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(7))
    
    # If on minor diagonal, must match existing diagonal letters
    if row + col == 6:
        diag_letters = set(grid[i][6-i] for i in range(7) if grid[i][6-i])
        if diag_letters:
            return list(diag_letters)
    
    valid = set('abcdefg') - used_row - used_col - {''}
    return list(valid)

def solve(grid, initial):
    empty = find_empty(grid)
    if not empty:
        return matches_initial(grid, initial)
    
    row, col = empty
    valid_letters = get_valid_letters(grid, row, col, initial)
    
    for letter in valid_letters:
        grid[row][col] = letter
        if (is_valid_row(grid, row) and 
            is_valid_col(grid, col) and 
            is_valid_minor_diagonal(grid)):
            if solve(grid, initial):
                return True
        grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['b','','','c','g','f','a'],
    ['','e','c','g','f','a','b'],
    ['e','c','g','f','a','b','d'],
    ['c','','f','','b','','e'],
    ['g','f','','b','d','',''],
    ['','','b','','e','c','g'],
    ['','','d','e','','g','f']
]

# Create working grid
grid = [row[:] for row in initial]

# Solve the puzzle
if solve(grid, initial):
    for row in grid:
        print(','.join(row))