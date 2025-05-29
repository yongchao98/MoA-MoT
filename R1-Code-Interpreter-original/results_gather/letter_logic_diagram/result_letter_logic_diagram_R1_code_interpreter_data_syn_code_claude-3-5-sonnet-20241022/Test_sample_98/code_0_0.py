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
        # Find any filled cell on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If cell is on minor diagonal, we need to match existing diagonal letter
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letters = grid[i][j]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['g','c','b','f','e','d',''],
    ['','','','e','d','a',''],
    ['','','e','d','a','g',''],
    ['f','e','','','g','c','b'],
    ['e','d','a','','','','f'],
    ['d','','g','c','','','e'],
    ['','g','','b','f','','']
]

if solve(grid):
    result = ""
    for row in grid:
        result += ','.join(row) + '\n'
    print(result.strip())
else:
    print("No solution exists")