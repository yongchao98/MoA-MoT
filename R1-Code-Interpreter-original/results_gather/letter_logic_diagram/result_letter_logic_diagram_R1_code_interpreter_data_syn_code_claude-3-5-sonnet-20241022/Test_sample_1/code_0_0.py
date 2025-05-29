def is_valid_row(grid, row):
    filled = [x for x in grid[row] if x]
    return len(filled) == len(set(filled))

def is_valid_col(grid, col):
    filled = [grid[r][col] for r in range(7) if grid[r][col]]
    return len(filled) == len(set(filled))

def is_valid_minor_diagonal(grid):
    diag = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
    return len(set(diag)) <= 1

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def is_valid_move(grid, row, col, letter):
    # Try the move
    grid[row][col] = letter
    
    # Check if valid
    valid = (is_valid_row(grid, row) and 
            is_valid_col(grid, col) and 
            is_valid_minor_diagonal(grid))
    
    # Undo the move
    grid[row][col] = ''
    return valid

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, must match existing diagonal letters
    if row + col == 6:
        diag_letters = set(grid[i][6-i] for i in range(7) if grid[i][6-i])
        if diag_letters:
            letters = diag_letters.pop()
    
    for letter in letters:
        if is_valid_move(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
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

# Solve the puzzle
if solve(initial):
    for row in initial:
        print(','.join(row))