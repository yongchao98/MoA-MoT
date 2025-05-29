def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if it's on minor diagonal
    if row + col == 6:  # Cell is on minor diagonal
        # Check if it matches other filled cells on minor diagonal
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

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, use the letter from another filled minor diagonal cell if exists
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letters = grid[i][j]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

if solve_grid(initial_grid):
    result = []
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")