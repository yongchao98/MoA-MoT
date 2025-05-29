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
    # Check if 'num' is valid in current position
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
    # If on minor diagonal, use the same letter as others on diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i]:
                letters = [grid[i][6-i]]
                break
        else:
            letters = 'abcdefg'
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    return False

# Initial grid
initial = [
    ['c','b','','','d','',''],
    ['b','g','f','','e','','c'],
    ['g','','d','e','a','c','b'],
    ['f','d','e','','c','b','g'],
    ['d','e','','','b','g',''],
    ['e','a','c','b','g','f','d'],
    ['','','b','g','','','e']
]

# Solve the puzzle
if solve(initial):
    for row in initial:
        print(','.join(row))
else:
    print("No solution exists")