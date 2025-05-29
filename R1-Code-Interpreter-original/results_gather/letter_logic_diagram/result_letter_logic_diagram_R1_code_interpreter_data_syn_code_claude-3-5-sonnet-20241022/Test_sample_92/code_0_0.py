def print_grid(grid):
    for row in grid:
        print(','.join(cell if cell else '' for cell in row))

def is_valid(grid, row, col, letter, size=7):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(size) if grid[i][col]]:
        return False
    
    # Check minor diagonal
    if row + col == size - 1:
        # If there's already a different letter on minor diagonal
        for i in range(size):
            if grid[i][size-1-i] and grid[i][size-1-i] != letter:
                return False
    
    return True

def find_empty(grid, size=7):
    for i in range(size):
        for j in range(size):
            if not grid[i][j]:
                return i, j
    return None

def solve_grid(grid, size=7):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = set('abcdefg')
    
    # If on minor diagonal, we need to check existing diagonal values
    if row + col == size - 1:
        for i in range(size):
            if grid[i][size-1-i]:
                letters = {grid[i][size-1-i]}
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
    ['a', 'f', 'b', 'd', '', 'c', ''],
    ['', '', 'd', '', 'c', '', 'a'],
    ['', '', 'e', '', '', '', 'f'],
    ['', '', 'c', '', '', 'f', ''],
    ['', '', '', 'a', '', 'b', ''],
    ['', '', '', 'f', 'b', '', ''],
    ['', 'a', 'f', '', '', 'e', '']
]

if solve_grid(initial_grid):
    print("<<<")
    print_grid(initial_grid)
    print(">>>")
else:
    print("No solution exists")