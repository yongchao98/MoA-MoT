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
    if row + col == 6:  # Minor diagonal position
        # Check if it matches other filled minor diagonal positions
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
    
    # If it's on minor diagonal, we need to match existing diagonal letter
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

# Initial grid
initial_grid = [
    ['', '', 'b', '', '', 'd', 'a'],
    ['', '', 'e', '', 'd', '', 'g'],
    ['', 'e', 'c', 'd', '', '', ''],
    ['', 'c', '', '', '', '', ''],
    ['', 'd', 'a', '', 'f', '', 'e'],
    ['', '', 'g', '', '', '', ''],
    ['', '', 'f', '', 'e', 'c', 'd']
]

if solve_grid(initial_grid):
    result = []
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")