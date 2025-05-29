def is_valid_move(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal - all must be the same letter
    if row + col == 6:  # Cell is on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal cells
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return i, j
    # Then fill other cells
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
    
    # If on minor diagonal, we need to match existing diagonal letter
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter:
            letters = diagonal_letter
    
    for letter in letters:
        if is_valid_move(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['a','e','f','g','','c','d'],
    ['e','f','g','b','','d','a'],
    ['','g','b','','','','e'],
    ['','b','','','a','','f'],
    ['','c','d','','e','','g'],
    ['c','d','a','e','f','g','b'],
    ['d','a','e','f','','b','']
]

if solve_grid(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")