def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal constraint - must be 'd' (from initial grid)
    if row + col == 6:  # If cell is on minor diagonal
        if letter != 'd':
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal if not complete
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then fill other cells
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
    
    # If on minor diagonal, only try 'd'
    if row + col == 6:
        letters = ['d']
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', 'a', 'c', '', '', 'e', 'd'],
    ['a', '', 'f', 'g', '', 'd', 'b'],
    ['c', 'f', '', '', '', '', ''],
    ['', 'g', 'e', 'd', '', '', ''],
    ['', 'e', 'd', 'b', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['d', '', 'a', '', '', '', '']
]

if solve(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")