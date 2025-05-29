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
    
    # If cell is on minor diagonal, we need to use the same letter as other filled cells
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letters = grid[i][j]
                break
        if letters == 'abcdefg':  # If no filled cell found on minor diagonal
            letters = 'a'  # We can use any letter
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'd', '', '', 'b', 'e', 'g'],
    ['', '', '', 'b', '', 'g', ''],
    ['f', '', '', 'e', 'g', '', ''],
    ['a', '', '', '', '', 'd', ''],
    ['', '', 'g', '', '', '', 'a'],
    ['', 'g', '', 'd', '', '', ''],
    ['', 'c', 'd', '', 'a', 'b', '']
]

if solve(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")