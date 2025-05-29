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
        # Check all minor diagonal cells that are filled
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
    
    # If cell is on minor diagonal, we must use the same letter as other filled cells on diagonal
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

# Initialize grid
initial_grid = [
    ['', 'a', 'e', 'g', 'f', 'b', ''],
    ['a', '', 'g', 'f', '', '', ''],
    ['', '', '', '', 'c', 'd', 'a'],
    ['', 'f', 'b', 'c', '', 'a', ''],
    ['', 'b', 'c', '', 'a', '', 'g'],
    ['b', '', '', '', '', 'g', ''],
    ['c', '', '', 'e', 'g', '', '']
]

# Solve the puzzle
if solve(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")