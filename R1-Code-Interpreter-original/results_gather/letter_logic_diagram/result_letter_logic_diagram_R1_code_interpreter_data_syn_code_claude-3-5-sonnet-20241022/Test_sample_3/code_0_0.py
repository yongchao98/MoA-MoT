def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this cell is part of it
    if row + col == 6:  # Cell is on minor diagonal
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

def solve_grid(grid):
    empty = find_empty(grid)
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
            if solve_grid(grid):
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

# Solve the puzzle
if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")